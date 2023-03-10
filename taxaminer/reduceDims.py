#!/usr/bin/env python

"""compute PCA of gene desprictors



Expects path to config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

import numpy as np

from . import classes
from . import checkInput

# import pathlib # to create directories
# import operator # for quick comparisons
# import scipy.stats as stats # for Pearson's R
# from itertools import product as itertools_product # to generate all possible oligonucleotides from base alphabet
# from Bio.Seq import Seq as BioPython_Seq # to count oligonucleotides (also overlapping ones! not implemented in normal count)

import numpy as np
import sys  # parse command line arguments
import logging
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import math
import umap

# prevents message "Loading [MathJax]/extensions/MathMenu.js" showing in plots
import plotly
plotly.io.kaleido.scope.mathjax = None


def distance(p1, p2):
    return math.sqrt(sum((p2[i] - p1[i]) ** 2 for i in range(len(p1))))


def farthest(pt, others):
    return max(others, key=lambda i: distance(pt, i))


def create_scree_plot(cfg, pca, plot_template):
    exp_var_cumul = np.cumsum(pca.explained_variance_ratio_)

    scree_plot = go.Figure(data=go.Scatter(
        x=[f'PC {i}' for i in range(1, pca.n_components_ + 1)],
        y=exp_var_cumul * 100))

    scree_plot.update_layout(
        title="Cummulative explained variance",
        xaxis_title="Principal Component",
        yaxis_title="explained variance [%]",
        template=plot_template,
    )
    scree_plot.write_image(f"{cfg.output_path}PCA/scree_plot.pdf")


def create_heatmap(cfg, pca, variables, plot_template):
    x = variables
    # TODO: adjust how many PCs to plot: use [:3, :]
    z = np.array(pca.components_)
    y = [f'PC {i}' for i in range(1, len(z) + 1)]

    # TODO: (custom) order of PCA variables
    order = variables

    heatmap = go.Figure(data=go.Heatmap(
        z=z,
        x=x,
        y=y,
        hoverongaps=False,
        colorscale='RdBu',
        reversescale=True,
    )

    )
    heatmap.update_layout(
        title="Contribution of variables",
        xaxis_title="variables",
        yaxis_title="Principal Component",
        xaxis_categoryarray=order,
        legend_title="contribution",
        template=plot_template
    )
    heatmap.update_xaxes(
        tickangle=-45,
    )

    heatmap.write_image(f"{cfg.output_path}PCA/heatmap.pdf")


def create_bar_plot(cfg, pca, plot_template):
    bar_plot = go.Figure(data=go.Bar(
        x=[f'PC {i}' for i in range(1, pca.n_components_ + 1)],
        y=pca.explained_variance_ratio_ * 100,
    )
    )
    bar_plot.update_layout(
        title="Explained variance",
        xaxis_title="Principal Component",
        yaxis_title="explained variance [%]",
        template=plot_template
    )

    bar_plot.write_image(f"{cfg.output_path}PCA/bar_plot.pdf")


def create_2d_biplot(cfg, pca, components, features, plot_template):
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

    # scaling of the loading to the data
    coordinates_as_tuples = [(row[0], row[1], row[2]) for row in components[:, :3]]
    farthest_point = farthest((0, 0, 0), coordinates_as_tuples)
    factor = max([np.abs(farthest_point[i]) / max(np.abs(loadings[:, i])) for i in range(len(farthest_point))])

    biplot = go.Figure()
    biplot.add_trace(go.Scatter(
        x=components[:, 0], y=components[:, 1], mode='markers', opacity=0.25, showlegend=False
    ))

    # random placement of the annotations to enhance readability
    xanchor = ["left", "center", "right"]
    yanchor = ["top", "middle", "bottom"]

    for i, feature in enumerate(features):
        biplot.add_shape(
            type='line',
            x0=0, y0=0,
            x1=loadings[i, 0] * factor,
            y1=loadings[i, 1] * factor,
        )
        biplot.add_annotation(
            x=loadings[i, 0] * factor,
            y=loadings[i, 1] * factor,
            ax=0, ay=0,
            xanchor=np.random.choice(xanchor),
            yanchor=np.random.choice(yanchor),
            showarrow=True,
            text=feature,
            yshift=5,
            xshift=5,
            arrowhead=1,
        )
    biplot.add_trace(go.Scatter(
        x=loadings[:, 0] * factor,
        y=loadings[:, 1] * factor,
        mode='markers',
        name='',
        hovertext=features,
        hoverinfo='text',
        marker=dict(
            size=2,
        ),
        showlegend=False
    )
    )

    biplot.update_layout(
        title="Biplot of PCA coordinates and contribution of variables",
        xaxis_title="PC 1",
        yaxis_title="PC 2",
        template=plot_template
    )

    biplot.write_image(f"{cfg.output_path}PCA/biplot.pdf")


def create_3d_biplot(cfg, pca, components, features, plot_template):
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

    # scaling of the loading to the data
    coordinates_as_tuples = [(row[0], row[1], row[2])
                             for row in components[:, :3]]
    farthest_point = farthest((0, 0, 0), coordinates_as_tuples)
    factor = max([np.abs(farthest_point[i]) / max(np.abs(loadings[:, i]))
                  for i in range(len(farthest_point))])

    biplot = go.Figure()
    biplot.add_trace(go.Scatter3d(
        x=components[:, 0], y=components[:, 1], z=components[:, 2],
        mode='markers', marker=dict(size=5, opacity=.25),
        name = "genes"
    ))

    coordinates_for_lines = pd.DataFrame(
        [list(loadings[int(i / 2)][:3]) if i % 2 == 0 else [0, 0, 0]
         for i in range(0, (len(loadings) * 2) - 1)],
        index=[features[int(i / 2)] if i % 2 == 0 else None
               for i in range(0, (len(features) * 2) - 1)])

    size_for_lines = pd.DataFrame([5 if i % 2 == 0 else 0
                                   for i in range(0, (len(features) * 2) - 1)])

    biplot.add_trace(go.Scatter3d(
        x=coordinates_for_lines[0] * factor,
        y=coordinates_for_lines[1] * factor,
        z=coordinates_for_lines[2] * factor,
        line=dict(
            color='black',
            width=2
        ),
        marker=dict(
            size=size_for_lines,
            color='black',
            # colorscale='Viridis',
        ),
        name='contribution',
        hovertext=coordinates_for_lines.index,
        hoverinfo='text',
    ))

    biplot.update_layout(
        title="Biplot of PCA coordinates and contribution of variables",
        scene=dict(
            xaxis_title="PC 1",
            yaxis_title="PC 2",
            zaxis_title="PC 2"
        ),
        template=plot_template
    )

    biplot.write_html(f"{cfg.output_path}PCA/biplot3d.html")

def write_contrib_report(cfg, pca, variables):

    contribution = pd.DataFrame(pca.components_.transpose(), index=variables)

    contribution.to_csv(f"{cfg.output_path}PCA/contribution_of_variables.csv",
                        header=[f'PC {i}' for i in range(1, pca.n_components_ + 1)])

def write_pca_summary(cfg, pca, variables):

    pca_summary = pd.DataFrame([[round(x, 4)
                                  for x in pca.explained_variance_ratio_],
                                [round(x, 4)
                                  for x in np.cumsum(pca.explained_variance_ratio_)]],
                               index=["Explained variance",
                                      "Cummulative explained variance"])

    pca_summary.to_csv(f"{cfg.output_path}PCA/pca_summary.csv",
                         header=[f'PC {i}'
                                 for i in range(1, pca.n_components_ + 1)])



def compute_umap(data):
    reducer = umap.UMAP(n_components=3)
    embedding = reducer.fit_transform(data)

    return embedding



def compute_pca(cfg):
    gene_info = pd.read_csv(f'{cfg.output_path}gene_info/imputed_gene_table.csv',
                            header=0, index_col=0)

    # variables to be used for PCA specified by user in config file
    pca_variables = cfg.input_variables.split(',')
    # select only the subset of variables
    data_w_nans = gene_info[pca_variables]

    # check if only one contig annotated
    if data_w_nans['c_name'].nunique() == 1:
        # drop the column related descriptors (no variance contained)
        contig_columns = ["c_num_of_genes",
                          "c_len",
                          "c_pct_assemby_len",
                          "c_pearson_r",
                          "c_pearson_p",
                          "c_gc_cont",
                          "c_gcdev",
                          "c_genelenm",
                          "c_genelensd"]
        contig_cov_cols = ["c_cov",
                           "c_covsd",
                           "c_covdev",
                           "c_genecovm",
                           "c_genecovsd"]
        indexed_cov_cols = [f'{col}_{index}' for col in contig_cov_cols for index in cfg.bam_paths.keys()]

        with open(f'{cfg.output_path}PCA/variables_excluded_from_PCA.txt',
                  'w') as excluded_cols_file:
            excluded_cols_file.write('contig related variables excluded due to single contig assembly:\n')
            excluded_cols_file.write(
                '\n'.join([col for col in (contig_columns + indexed_cov_cols) if col in data_w_nans.columns]))
            excluded_cols_file.write('\n')

        data = data_w_nans.drop(columns=contig_columns + indexed_cov_cols, errors='ignore')
    else:
        data = data_w_nans

    # CLEAN DATA (FROM NaNS) -- COLUMNWISE
    # if there are any columns with at least 30% NaN values, save them to a new data frame
    nan_col_count = data.isna().sum()
    row_count = len(data.index)
    nan_columns = pd.DataFrame(index=data.index)
    for col, count in nan_col_count.items():
        if (count / row_count) > .3:
            nan_columns.insert(0, 'new_col', data[col])
            nan_columns.rename(columns={"new_col": col}, inplace=True)
            # keep only those columns with less than 30% NaN values
            data = data.drop(columns=col)
    with open(f'{cfg.output_path}PCA/variables_excluded_from_PCA.txt', 'a') as excluded_cols_file:
        excluded_cols_file.write('contig related variables excluded due to more then 30% NaNs:\n')
        excluded_cols_file.write('\n'.join(list(nan_columns.columns)))
        excluded_cols_file.write('\n')

    # CLEAN DATA FROM NaNS -- ROWWISE
    # save the genes containing NaN values to a specific data frame
    genes_w_nans = data[data.isna().any(axis=1)]
    # and log their IDs
    if not genes_w_nans.empty:
        logging.info(f'Following genes with null values are excluded from PCA:\n{genes_w_nans.index.values}')

    # keep working only with the genes without NaNs (complete rows / cases)
    data = data.dropna()

    # TODO: filter coverage modes

    # :::::::::::::::::::: PCA ::::::::::::::::::::::::::::::

    scaled_data = StandardScaler().fit_transform(data.iloc[:, 1:])
    pca = PCA()
    components = pca.fit_transform(scaled_data)

    pca_coordinates = pd.DataFrame(data=components[:, :3],
                                   columns=['PC 1', 'PC 2', 'PC 3'],
                                   index=data.index)

    # umap_embedding = compute_umap(scaled_data)
    # umap_df = pd.DataFrame(data=umap_embedding
    #                     , columns=['PC 1', 'PC 2', 'PC 3'],index=data.index)

    plot_template = dict(
        layout=go.Layout(
            template="plotly_white",
            font=dict(
                size=14
            ),
            xaxis=dict(automargin=True, title_standoff=15),
            yaxis=dict(automargin=True, title_standoff=5)
        ))

    create_2d_biplot(cfg, pca, components, data.columns[1:], plot_template)
    create_heatmap(cfg, pca, data.columns[1:], plot_template)
    create_scree_plot(cfg, pca, plot_template)
    create_bar_plot(cfg, pca, plot_template)
    create_3d_biplot(cfg, pca, components, data.columns[1:], plot_template)

    write_contrib_report(cfg, pca, data.columns[1:])
    write_pca_summary(cfg, pca, data.columns[1:])

    return pca, pca_coordinates, data.columns[1:]

    # TODO: horns parallel analysis


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    compute_pca(cfg)


if __name__ == '__main__':
    main()
