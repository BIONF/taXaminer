#!/usr/bin/env python

"""Create and modifiy output

Modifies the HTML 3D plot file to make is self-contained and
the text in the hoverwindow selectable. Creates final 3D and Krona plot

Expects processed config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"




from bs4 import BeautifulSoup as bs
from jsmin import jsmin

import json

from . import checkInput
from . import reduceDims
from . import compFeatures

import taxopy
import sys
import csv
import pathlib
import subprocess
import logging
import pandas as pd
from Bio import SeqIO

import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.colors import n_colors
import plotly.io as pio

# mathjax is not required for saving of the plots;
# it does though leave a message bar on the plots
pio.kaleido.scope.mathjax = None

def rotate_z(x, y, z, theta):
    """
    Calculate 3d coordinates for a circular rotation of the scatterplot
    :param x: x coordinate
    :param y: y coordinate
    :param z: z coordinate
    :param theta: rotation step
    :return: triple of x,y,z coordinates
    """
    c = x + 1j * y
    return np.real(np.exp(1j * theta) * c), np.imag(np.exp(1j * theta) * c), z


def gene_data2panda(cfg, assignment_df, pca_coordinates):
    gene_features = pd.read_csv(f"{cfg.output_path}gene_info/raw_gene_table.csv",
                                index_col=0)

    df = gene_features.merge(
        pca_coordinates, left_index=True, right_index=True).merge(
        assignment_df, left_index=True, right_index=True)

    return df


def insert_string2string(string, insert_string, interval):
    new_string = ""
    for i in range(-(-len(string) // interval)):
        new_string = new_string + string[i * interval:interval * (
                    i + 1)] + insert_string
    return new_string


def assess_superkingdom(plot_df):
    # get unique plot labels
    plot_labelIDs = plot_df['plot_labelID'].unique()
    # identify their respective superkingdom
    for id in plot_labelIDs:
        if not pd.isna(id):
            taxon = taxopy.Taxon(id, TAX_DB)
            superkingdom = taxon.rank_name_dictionary.get('superkingdom')
            plot_df.loc[
                plot_df['plot_labelID'] == id, 'superkingdom'] = superkingdom
        else:
            plot_df.loc[plot_df['plot_labelID'] == id, 'superkingdom'] = np.nan


def colourscale_to_lengths(lengths, colourscale):
    if min(lengths[0]) == max(lengths[0]):
        normalized_lengths = [1 for x in lengths[0]]
    else:
        normalized_lengths = [
            compFeatures.rescale(x, min(lengths[0]), max(lengths[0]), 0.3,
                                      1, np.mean(lengths[0])) for x in
            lengths[0]]

    colouring = pd.DataFrame(
        px.colors.sample_colorscale(colourscale, normalized_lengths),
        index=lengths.index)
    return colouring


def default_colouring(plot_df, query_label, sk_sorted_labels, colourscale):

    # non-taxonomy informed colorscale
    distinct_labels_count = plot_df['plot_label'].value_counts(
        dropna=False).to_frame()
    distinct_labels_count.drop(['Unassigned', query_label], inplace=True,
                               errors='ignore')
    # distinct_labels_count = pd.DataFrame(
    #     {'plot_label_count': distinct_labels_count['plot_label']},
    #     index=distinct_labels_count.index)

    # sk_df = pd.DataFrame(
    #     plot_df[["superkingdom", "plot_label"]]
    # ).set_index("plot_label")
    # sk_df = sk_df[~sk_df.index.duplicated(keep='first')]
    # distinct_labels_count = distinct_labels_count.join(sk_df)
    # distinct_labels_count = distinct_labels_count.sort_values(by='superkingdom')

    sort_index = dict(zip(sk_sorted_labels, range(len(sk_sorted_labels))))
    distinct_labels_count['sk_sort'] = distinct_labels_count.index.map(sort_index)
    distinct_labels_count.sort_values('sk_sort', inplace=True)
    distinct_labels_count.drop('sk_sort', axis=1, inplace=True)

    if not distinct_labels_count.empty:
        marker_style = pd.DataFrame(
            {"plot_colour": px.colors.sample_colorscale(
                colourscale, np.arange(0, 1, 1 / len(distinct_labels_count))),
             "plot_opacity": [1 for i in range(distinct_labels_count.shape[0])]},
            index=distinct_labels_count.index)
        plot_df = plot_df.join(marker_style, on='plot_label')

    plot_df.loc[plot_df[
                    'plot_label'] == query_label, 'plot_colour'] = 'rgb(159,161,179)'
    plot_df.loc[plot_df[
                    'plot_label'] == "Unassigned", 'plot_colour'] = 'rgb(211,213,220)'

    plot_df.loc[plot_df['plot_label'] == query_label, 'plot_opacity'] = 0.8
    plot_df.loc[plot_df['plot_label'] == "Unassigned", 'plot_opacity'] = 0.65

    return plot_df


def labels_sorted_sk(plot_df):

    assess_superkingdom(plot_df)

    sk_df = pd.DataFrame(
        plot_df[["superkingdom", "plot_label"]]
    ).set_index("plot_label")
    sk_df.index.name = 'plot_label'
    sk_df = sk_df[~sk_df.index.duplicated(keep='first')]
    sk_df = sk_df.sort_values(by=['superkingdom', 'plot_label'])

    return list(sk_df.index)



def save_gif(cfg, fig):
    # save a gif of the 3D plot
    new_fig = fig.update_layout(showlegend=False,
                            font=dict(size=10, color="grey"))
    gif = GIF()
    three_d_scatter_rotate(gif, new_fig, 60,
                            gif_kwargs={'length': 45000})
    # custom path in gif_kwargs does not work due to bug
    # 'gif_path': f'{cfg.output_path}taxonomic_assignment/3D_plot.gif'



def create_3D_plot(cfg, all_data_df, pca_obj, variables, pca_coordinates, query_label):
    """

    Args:
      cfg:
      pca_coordinates:

    Returns:

    """



    # ________________ PLOT PREPARATION _______________ #

    plot_df = all_data_df

    # frequency information for labels in plots
    plot_df['label_count'] = plot_df['plot_label'].map(
        plot_df['plot_label'].value_counts())
    # append frequency information to label
    plot_df['plot_label_freq'] = plot_df['plot_label'] + ' (' + plot_df[
        'label_count'].astype(str) + ')'

    # list of the plot labels, sorted by superkindgom
    # and alphabetically within superkingdom
    sk_sorted_labels = labels_sorted_sk(plot_df)

    ## add sequence information
    # read proteins fasta file
    # delete everything after first whitespace in fasta header to match table
    # add html line breaks to sequence so hover windows stay narrrow
    fasta_sequences = {rec.id.split()[0]: rec.seq for rec in
                       SeqIO.parse(cfg.proteins_path, "fasta")}
    # merge to data table
    plot_df['prot_seq'] = plot_df['fasta_header'].map(fasta_sequences)
    # reformat prot seq for plot hover label (insert html line breaks)
    plot_df['plot_seqs'] = plot_df['prot_seq'].apply(
        lambda x: insert_string2string(str(x), '</br>', 70))

    ## prepare data for hover window display
    # replace bools to yes and no
    plot_df.loc[plot_df['g_terminal'] == 0, 'g_terminal'] = 'no'
    plot_df.loc[plot_df['g_terminal'] == 1, 'g_terminal'] = 'yes'
    # all gene coverage information in one column
    g_cov_cols = plot_df.filter(regex=("g_cov_[0-9]*"))
    plot_df['g_coverages'] = g_cov_cols.apply(
        lambda x: '; '.join(x.dropna().astype(str)), axis=1)
    # same for gene coverage deviation
    g_covdev_cols = plot_df.filter(regex=("g_covdev_c_[0-9]*"))
    plot_df['g_covdeviations'] = g_covdev_cols.apply(
        lambda x: '; '.join(x.dropna().astype(str)), axis=1)




    traces = plot_df['plot_label'].unique()
    traces_reordered = [plot_df.loc[plot_df['plot_label'] == query_label,
                            'plot_label'][0]]
    try:
        traces_reordered.append(plot_df.loc[plot_df['plot_label'] == "Unassigned",
                                    'plot_label'][0])
    except:
        pass

    remaining_traces = np.setdiff1d(traces, traces_reordered)
    sorted_remaining_traces = [trace for trace in sk_sorted_labels if
                               trace in remaining_traces]
    traces_reordered.extend(sorted_remaining_traces)

    sort_index = dict(zip(traces_reordered, range(len(traces_reordered))))
    plot_df['sk_sort'] = plot_df['plot_label'].map(sort_index)
    plot_df.index.name = "g_name"
    plot_df.sort_values(['sk_sort'], inplace=True)
    plot_df.drop('sk_sort', axis=1, inplace=True)

    ## assign colours
    plot_df = default_colouring(plot_df, query_label, sk_sorted_labels, 'Rainbow')

    ## create plot
    # subset data into three groups to load them in distinct traces
    # this enables to stack the data points in the desired order in the plot
    # i.e. putting the less interesting (background) data points with query assignment in the background

    fig = go.Figure()
    for label in traces_reordered:
        labels = plot_df.loc[plot_df['plot_label'] == label, :]
        fig.add_trace(go.Scatter3d(
            x=plot_df.loc[plot_df['plot_label'] == label, 'PC 1'],
            y=plot_df.loc[plot_df['plot_label'] == label, 'PC 2'],
            z=plot_df.loc[plot_df['plot_label'] == label, 'PC 3'],
            name=plot_df.loc[plot_df['plot_label'] == label, 'plot_label_freq'][0],
            mode='markers',
            marker=dict(
                size=2,
                color=
                plot_df.loc[plot_df['plot_label'] == label, 'plot_colour'][0],
                opacity=plot_df.loc[
                    plot_df['plot_label'] == label, 'plot_opacity'][0],
            ),
            text=[f"</br>ID: {index} | Scaffold: {row['c_name']} \
                            </br>Taxonomic assignment: {row['taxon_assignment']} | Label: {row['plot_label']} \
                            </br>Coverage: {row['g_coverages']} (SD from contig mean: {row['g_covdeviations']}) \
                            </br>Terminal: {row['g_terminal']} (Genes on contig: {row['c_num_of_genes']}) \
                            </br>LCA: {row['lca']} \
                            </br>Best hit: {row['best_hit']} (e-value: {row['bh_evalue']}; pident: {row['bh_pident']}; bitscore: {row['bh_bitscore']}) \
                            </br>Seq: {row['plot_seqs']}" for index, row in plot_df.loc[plot_df['plot_label'] == label].iterrows()],
            hoverinfo='text',
            showlegend=True
            # legendgroup="genes",
            # legendgrouptitle_text="Taxonomic assignment"
        ))

    # tight layout
    fig.update_layout(
        template="plotly_white",
        title=cfg.taxon_name,
        scene=dict(
            xaxis_title="PC 1",
            yaxis_title="PC 2",
            zaxis_title="PC 3"
        ),
        scene_camera_eye=dict(x=-1.9, y=1.7, z=1.45),
        font=dict(
            size=14
        ),
        legend_title_text='Taxonomic assignment',
        margin=dict(
            l=10,
            r=10,
            b=10,
            t=40
        ),
        legend={'itemsizing': 'constant',
                'groupclick': "toggleitem",
                'tracegroupgap': 2}
    )



    # pio.kaleido.scope.default_width = 700
    # pio.kaleido.scope.default_height = 500
    fig.write_image(cfg.output_path + "taxonomic_assignment/3D_plot.pdf",
                    engine='kaleido',
                    width = 850,
                    height = 650)


    #save_gif(cfg,fig)

    ### ADD PCA CONTRIBUTION OF VARIABLES ###
    loadings = pca_obj.components_.T * np.sqrt(pca_obj.explained_variance_)
    # scaling of the loading to the data
    coordinates_as_tuples = [(row[1]['PC 1'], row[1]['PC 2'], row[1]['PC 3'])
                             for row in pca_coordinates.iterrows()]
    farthest_point = reduceDims.farthest((0, 0, 0), coordinates_as_tuples)
    factor = max(
        [np.abs(farthest_point[i]) / max(np.abs(loadings[:, i])) for i in
         range(len(farthest_point))]) * 0.5

    coordinates_contribution = pd.DataFrame({
        'center': [(0, 0, 0) for i in range(loadings.shape[0])],
        'loading': [(row[0], row[1], row[2]) for row in loadings]},
        index=variables)

    lengths = pd.DataFrame(
        [reduceDims.distance(row[1]['center'], row[1]['loading']) for row in
         coordinates_contribution.iterrows()],
        index=variables)
    colour_lengths = colourscale_to_lengths(lengths, 'Greys')

    coordinates_contribution = coordinates_contribution.join(colour_lengths)

    for variable in coordinates_contribution.iterrows():
        x = [variable[1]['center'][0], variable[1]['loading'][0] * factor]
        y = [variable[1]['center'][1], variable[1]['loading'][1] * factor]
        z = [variable[1]['center'][2], variable[1]['loading'][2] * factor]
        fig.add_trace(go.Scatter3d(
            x=x,
            y=y,
            z=z,
            visible="legendonly",
            legendgroup="variables",
            legendgrouptitle_text="Contribution<br>of variables",
            line=dict(
                color=variable[1][0],
                width=2
            ),
            marker=dict(
                size=[0, 5],
                color=variable[1][0]
            ),
            name=variable[0],
            hovertext=variable[0],
            hoverinfo='text'
        ))

    fig.update_layout(legend=dict(groupclick="togglegroup"),
                      modebar_add=["v1hovermode"])
    fig.write_html(cfg.output_path + "taxonomic_assignment/3D_plot.html")

    return plot_df


def create_krona_plot(cfg, data):
    krona_df = data['taxon_assignmentID']
    krona_df.to_csv(cfg.output_path + 'tmp/krona.tsv', sep='\t', index=True,
                    index_label='#g_name', header=['taxID'])
    cmd_krona = f'{cfg.krona} {cfg.output_path}/tmp/krona.tsv -o {cfg.output_path}taxonomic_assignment/krona.html'
    out_krona = subprocess.run([cmd_krona], shell=True, capture_output=True)
    if out_krona.returncode != 0:
        logging.error(f'creation of krona plot failed:\n{out_krona}')
        logging.error('Error message:\n' + out_krona.stderr.decode())


def create_gene_table_taxon_assignment(cfg, gff_df, df):

    out_df = df.merge(gff_df[['start', 'end', 'upstream_gene', 'downstream_gene']],
                      left_index=True, right_index=True)

    columns_renamed = {}
    for col in out_df.columns:
        if 'PC ' in col:
            columns_renamed[col] = f"Dim.{col.split()[1]}"
    out_df.rename(columns=columns_renamed, inplace=True)
    out_df.to_csv(f"{cfg.output_path}taxonomic_assignment/gene_table_taxon_assignment.csv", index_label='g_name')


def create_plots(cfg, genes, pca_obj, variables, pca_coordinates,
                 query_label, tax_db, gff_df):
    global TAX_DB
    TAX_DB = tax_db

    all_data_df = gene_data2panda(cfg, genes, pca_coordinates)
    create_gene_table_taxon_assignment(cfg, gff_df, all_data_df)
    create_3D_plot(cfg, all_data_df, pca_obj, variables, pca_coordinates, query_label)
    #TODO: krona optional
    # if cfg.create_krona:
    create_krona_plot(cfg, all_data_df)

    return all_data_df


def make_selectable(html_path):
    """Make text in hoverwindows in 3D plot selectable.

    Replaces parts of javascript code to make the text in the
    hoverwindow selectable. Highly specified function and may not work
    with other versions of htmlwidgets and plotly (developed for
    htmlwidgets-1.5.7 and plotly-4.9.0; compatiblity tested with
    htmlwidgets-1.5.4 and plotly-4.10.0)

    Args:
      html_path(str): path to javascript file where to replace code
    """

    with open(html_path, 'r') as file_in:
        filedata = file_in.read()

    old = ('r._paperdiv=r._container.selectAll(".svg-container").data([0]),'
        'r._paperdiv.enter().append("div").classed("user-select-none",!0).classed("svg-container",!0).style("position","relative"),'
        'r._glcontainer=r._paperdiv.selectAll(".gl-container").data([{}]),r._glcontainer.enter().append("div").classed("gl-container",!0),'
        'r._paperdiv.selectAll(".main-svg").remove(),r._paperdiv.select(".modebar-container").remove(),'
        'r._paper=r._paperdiv.insert("svg",":first-child").classed("main-svg",!0),r._toppaper=r._paperdiv.append("svg").classed("main-svg",!0),'
        'r._modebardiv=r._paperdiv.append("div"),delete r._modeBar,r._hoverpaper=r._paperdiv.append("svg").classed("main-svg",!0)')

    # check if text to be replaced is in current file
    if old not in filedata:
        # probable reason is other version of plotly/htmlwidgets
        logging.warning("Text in interactive plot could not be made selectable.\
                        Check versions of htmlwidgets and plotly")

    new = ('r._paperdiv=r._container.selectAll(".svg-container").data([0]),'
        'r._paperdiv.enter().append("div").classed("svg-container",!0).style("position","relative"),'
        'r._glcontainer=r._paperdiv.selectAll(".gl-container").data([{}]),r._glcontainer.enter().append("div").classed("gl-container",!0),'
        'r._paperdiv.selectAll(".main-svg").remove(),r._paperdiv.select(".modebar-container").remove(),'
        'r._paper=r._paperdiv.insert("svg",":first-child").classed("user-select-none",!0).classed("main-svg",!0),'
        'r._toppaper=r._paperdiv.append("svg").classed("user-select-none",!0).classed("main-svg",!0),'
        'r._modebardiv=r._paperdiv.append("div").classed("user-select-none",!0),delete r._modeBar,'
        'r._hoverpaper=r._paperdiv.append("svg").classed("user-select-none",!0).classed("main-svg",!0)')


    filedata = filedata.replace(old, new)

    with open(html_path, 'w') as file_out:
        file_out.write(filedata)


def make_selfcontained(html_file_path, html_dependencies_dir, html_out_path):
    """Make HTML file self-contained.

    Replaces links to scripts in the head part of HTML with
    script content to make it selfcontained

    Args:
      html_file_path(str): path to HTML plot file
      html_dependencies_dir(str): path to dir with HTML dependencies
      html_out_path(str): path to save self-contained HTML file to
    """

    with open(html_file_path) as html_file:
        soup = bs(html_file, 'html.parser')
        # find script elements in HTML file
        scripts = soup.find_all("script")
        for script_item in scripts:
            # check if element has src (link) attribute
            try:
                script_path = script_item.attrs['src']
            except:
                continue

            # script that is cruical to make text in hover selectable
            # and needs to be modified for that
            if 'plotly-main-' in script_path:
                make_selectable(html_dependencies_dir + script_path)

            # other scripts: read content and minify if it is javascript
            with open(html_dependencies_dir + script_path) as dependency_file:
                if script_path.endswith('.js'):
                    dependency_text = jsmin(dependency_file.read())
                else:
                    dependency_text = dependency_file.read()

            # remove script element from HTML
            script_item.extract()
            # add element with contend of script
            selfcontained_script = soup.new_tag('script')
            selfcontained_script.string = dependency_text
            soup.head.append(selfcontained_script)

        # replace links to css files with their content
        # elements are found by tag link
        head_links = soup.head.find_all('link')
        for link in head_links:
            # get path to css file
            css_path = link.attrs['href']

            # remove element from HTML
            link.extract()
            # retrieve content of css file
            with open(html_dependencies_dir + css_path) as css_file:
                styling_tag = soup.new_tag('style')
                styling_tag.string = css_file.read()

                # add content of css directly into HTML head
                soup.head.append(styling_tag)

    # write modified HTML to output location
    with open(html_out_path, 'w') as html_file:
        html_file.write(str(soup))


def change_title(cfg, output_path, html_out_path):
    """Change title of html file to species name.

    Reads species name from file temporarily created from
    taxonomic_assignment.py and replaces old title with species name

    Args:
      output_path: path to directory of taXaminer report
      html_out_path: path to HTML plot file
      :param cfg:
    """

    with open(html_out_path) as html_file:
        soup = bs(html_file, 'html.parser')
        # find and delete old title
        if soup.find("title"):
            title = soup.find("title")
            title.extract()
        # replace with new title
        new_title = soup.new_tag('title')
        new_title.string = cfg.taxon_name
        soup.head.append(new_title)

    # write modified HTML to same location
    with open(html_out_path, 'w') as html_file:
        html_file.write(str(soup))


def perform_adjustments(cfg):
    """Collectively perform all modifications to HTML plot file.

    HTML is made self-contained, text is made selectable and
    title is changed to species name. Location of HTML is changed from
    tmp directory to final location in taxonomic_assignment dir

    Args:
      cfg: Config object with config parameters
    """

    # temporary location of HTML
    html_file_path = cfg.output_path + 'taxonomic_assignment/3D_plot.html'
    html_dependencies_dir = cfg.output_path + 'tmp/'
    # final location
    html_out_path = cfg.output_path + 'taxonomic_assignment/3D_plot.html'

    # makes selfcontained and calls script to make text selectable
    # make_selfcontained(html_file_path, html_dependencies_dir, html_out_path)

    make_selectable(html_file_path)
    change_title(cfg, cfg.output_path, html_file_path)


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    #TODO: add genes
    genes = None
    perform_adjustments(cfg)
    create_plots(cfg, genes)


if __name__ == '__main__':
    main()
