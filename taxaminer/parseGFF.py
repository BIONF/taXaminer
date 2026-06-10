#!/usr/bin/env python

"""Prepare data for running taXaminer

Parses the GFF file

Prepares coverage data for further processing. Raw reads are mappend
and mapping files (BAM) are converted to per base coverage (PBC) files

Compute protein FASTA by parsing gene sequences baed on coordinates from
GFF and FASTA with scaffold sequences (assembly FASTA) and converting
it to amino acid sequence

Expects prepared config file with information regarding available
coverage information and corresponding paths to find the files
(preparation of config by prepare_and_check.py)
"""
__author__ = "Freya Arthen"

import logging
import sys
import pandas as pd
import copy


from . import checkInput


def check_gff(cfg):
    pass

def remove_prefix(text, prefix):
    """
    Removes given prefix from text, returns stripped text.

    Args:
      text(str): text to remove prefix from
      prefix(str): prefix to remove from text

    Returns:
      (str) stripped text
    """
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def strip_ID(id):
    """
    Removes attribute-linked prefixes from IDs, returns stripped ID.

    Args:
      id(str): ID with prefix specifying type of attribute

    Returns:
      (str) stripped ID

    """
    if pd.isna(id):
        return None

    prefix_types = ['gene', 'GENE', 'Gene',
                    'transcript', 'TRANSCRIPT', 'Transcript',
                    'rna', 'RNA', 'Rna',
                    'cds', 'CDS', 'Cds']
    prefix_connectors = [':', '-']
    prefix_list = [f"{type}{con}" for con in prefix_connectors for type in prefix_types]

    for prefix in prefix_list:
        if id.startswith(prefix):
            id = remove_prefix(id, prefix)
            return id
    return id


def attrs2dict(attr_list):
    """ Parse the GFF additional attributes list into a dict

    :param attr_list:
    :return:
    """

    return {attr.split('=')[0]: attr.split('=')[1] for attr in
            attr_list.split(';') if attr}


def spline2dict(spline, include_pseudogenes):

    spline_dict = {}
    # direct inference
    spline_dict['scaffold'] = spline[0]
    spline_dict['source'] = spline[1]
    if include_pseudogenes and spline[2] == 'pseudogene':
        spline_dict['type'] = 'gene'
    else:
        spline_dict['type'] = spline[2]
    spline_dict['start'] = int(spline[3])
    spline_dict['end'] = int(spline[4])
    spline_dict['score'] = spline[5]
    spline_dict['strand'] = spline[6]
    spline_dict['phase'] = spline[7]
    spline_dict['add_attrs'] = attrs2dict(spline[8])
    #spline_dict['add_attrs']['ID'] = strip_ID(spline_dict['add_attrs'].get('ID'))
    #spline_dict['add_attrs']['Parent'] = strip_ID(spline_dict['add_attrs'].get('Parent'))
    spline_dict['id'] = spline_dict['add_attrs'].get('ID')
    spline_dict['parent_id'] = spline_dict['add_attrs'].get('Parent')
    transl_table = spline_dict['add_attrs'].get('transl_table')
    spline_dict['transl_table'] = transl_table if transl_table else '1'

    # indirect inference
    spline_dict['gene_id'] = spline_dict['id'] if spline_dict['type'] == 'gene' else None
    spline_dict['length'] = spline_dict['end'] - spline_dict['start'] + 1
    if not spline_dict.get('id'):
        # no ID information only acceptable for features without subfeatures
        spline_dict['id'] = f"{spline_dict.get('parent_id')}-{spline_dict['type']}-{spline_dict['start']}"

    return spline_dict


def add_missing_ids(gff_df):
    # if features don't have ID attribute, concatenate type and line number
    no_id = gff_df['id'].isnull()
    gff_df.loc[no_id, 'id'] = gff_df['type'] + '-' + gff_df.index.astype(
        str)

    return gff_df


def remove_duplicate_ids(gff_df):
    # check for duplicate gene IDs
    if gff_df.loc[gff_df['type'] == 'gene'].duplicated(subset='id',
                                                       keep=False).any():
        duplicate_genes = gff_df.loc[gff_df['type'] == 'gene'].loc[
            gff_df.loc[gff_df['type'] == 'gene'].duplicated(
                subset='id', keep=False)]
        duplicate_ids = duplicate_genes['id'].unique()
        # TODO: include partial genes?
        # for dup_gene in duplicate_genes.itertuples():
        #     if 'part' in dup_gene.add_attrs.keys():

        logging.warning(f"Detected {len(duplicate_ids)} "
                        f"duplicate gene IDs. Will procede without them.\n"
                        f"Duplicate IDs: {duplicate_ids}")
        gff_df.drop(gff_df.loc[gff_df['id'].isin(duplicate_ids)].index,
                    inplace=True)
        gff_df.drop(gff_df.loc[gff_df['gene_id'].isin(duplicate_ids)].index,
                    inplace=True)

    # append line number to duplicate ids
    duplicate_index = gff_df.duplicated(subset='id', keep=False)
    gff_df.loc[(duplicate_index) & (gff_df['type'] != 'gene'), 'id'] = \
        gff_df['id'] + '-' + gff_df.index.astype(str)

    return gff_df

def set_ids2index(gff_df):
    gff_df.set_index('id', inplace=True)
    return gff_df


def strip_ids(gff_df):

    columns_w_ids = ['id', 'gene_id', 'transcript_id', 'parent_id', 'coding_features', 'unique_coding_id']

    for col in columns_w_ids:
        gff_df[col] = gff_df.apply(lambda row: [strip_ID(val) for val in row[col]] if type(row[col]) == list else strip_ID(row[col]), axis=1)

    return gff_df


def ident_parent_path(feature, search_features):
    parent_id = [feature.get('parent_id')]
    updated = True
    while updated:
        updated = False
        for feature in search_features:
            if feature.get('id') == parent_id[-1]:
                parent_id.append(feature.get('parent_id'))
                updated = True
                break

    return parent_id

def multi_transcript_gene(non_transcript_types, gene_features, transcripts):
    coding_features = None
    if 'exon' in non_transcript_types:
        if 'CDS' in non_transcript_types:
            # CDS is coding type
            all_coding_features = [cds for cds in gene_features if
                                   cds.get('type') == 'CDS']
        else:
            # exon is coding type
            all_coding_features = [exon for exon in gene_features if
                                   exon.get('type') == 'exon']
    elif 'CDS' in non_transcript_types:
        # CDS is coding type
        all_coding_features = [cds for cds in gene_features if
                               cds.get('type') == 'CDS']
    else:
        # mRNA is coding type
        coding_features = sorted(transcripts,
                                 key=lambda dict: dict.get('length'))[-1:]
        longest_transcript = coding_features[0]


    if not coding_features:
        if len(transcripts) == 1:
            longest_transcript = transcripts[0]
            coding_features = sorted(all_coding_features,
                                     key=lambda dict: dict.get(
                                         'start'))
        else:
            # identify longest transcript subfeatures
            gene_cds_len = 0
            unsorted_coding_features = None
            longest_transcript = None
            for transcript in transcripts:
                associated_coding_features = \
                    [feat for feat in all_coding_features
                     if transcript.get('id') in ident_parent_path(feat,
                                                                  transcripts + gene_features)]
                transcript_cds_len = sum(
                    [feat.get('length') for feat in associated_coding_features])
                if transcript_cds_len > gene_cds_len:
                    gene_cds_len = transcript_cds_len
                    unsorted_coding_features = associated_coding_features
                    longest_transcript = transcript
            if unsorted_coding_features:
                coding_features = sorted(unsorted_coding_features,
                                         key=lambda dict: dict.get('start'))
            else:
                return None, None

    return coding_features, longest_transcript

def process_gene(gene, parent_ids, transcripts, gene_features, non_transcript_types):


    gene_biotype = 'mRNA' # = default, will otherwise be replaced
    longest_transcript = None

    if transcripts:
        transcript_types = [transcript.get('type') for transcript in
                            transcripts]
        unique_transcript_types = set(transcript_types)
        if len(unique_transcript_types) == 1:
            transcript_type = transcript_types[0]

            if transcript_type == 'mRNA':
                coding_features, longest_transcript = multi_transcript_gene(
                    non_transcript_types, gene_features, transcripts)

            elif transcript_type == 'exon':
                if 'CDS' in non_transcript_types:
                    # CDS is coding type
                    unsorted_coding_features = [cds for cds in gene_features if cds.get('type') == 'CDS']
                    coding_features = sorted(unsorted_coding_features,
                                             key=lambda dict : dict.get('start'))
                else:
                    # exon is coding type
                    unsorted_coding_features = [exon for exon in transcripts]
                    coding_features = sorted(unsorted_coding_features,
                                                   key=lambda dict: dict.get(
                                                       'start'))

            elif transcript_type == 'CDS':
                # CDS is coding type
                unsorted_coding_features = [cds for cds in transcripts]
                coding_features = sorted(unsorted_coding_features,
                                               key=lambda dict: dict.get(
                                                   'start'))
            else:
                # non-protein-coding RNA
                gene_biotype = transcript_type
                gene['coding_features'] = None
                gene['coding_coordinates'] = None
                gene['transcript_id'] = None
                gene['gene_biotype'] = gene_biotype
                return [gene]

        elif unique_transcript_types == {'exon', 'CDS'}:
            # CDS is coding type
            unsorted_coding_features = [cds for cds in transcripts if
                                  cds.get('type') == 'CDS']
            coding_features = sorted(unsorted_coding_features,
                                     key=lambda dict: dict.get(
                                         'start'))
        elif 'mRNA' in unique_transcript_types:
            coding_transcripts = [feat for feat in transcripts if
                                  feat.get('type') == 'mRNA']
            if 'exon' in unique_transcript_types:
                non_transcript_types.append('exon')
                gene_features += [feat for feat in transcripts
                                  if feat.get('type') == 'exon']
            if 'CDS' in unique_transcript_types:
                non_transcript_types.append('CDS')
                gene_features += [feat for feat in transcripts if
                                  feat.get('type') == 'CDS']
            coding_features, longest_transcript = multi_transcript_gene(
                non_transcript_types, gene_features, coding_transcripts)

        else:
            logging.warning(f"Transcript type of gene {gene.get('id')} ambiguous."
                            f" Will add without coding features.")
            return [gene]


    else:
        # gene == coding type
        gene_cds = copy.deepcopy(gene)
        gene_cds['type'] = 'CDS'
        coding_features = [gene_cds]
        longest_transcript = gene_cds

    if not coding_features:
        logging.warning(f"Coding features of {gene.get('id')} inconclusive."
                        f" Will add without coding features.")
        return [gene]

    # identify unambigous id amongst coding features to use as potential
    # mapping value of protein headers
    unique_coding_id = None
    cf_ids = [feat.get('id') for feat in coding_features]
    # coding features all have same ID
    if len(set(cf_ids)) == 1:
        unique_coding_id = cf_ids[0]
        gene['unique_coding_id'] = unique_coding_id
        for i, feat in enumerate(coding_features):
            feat['id'] = f"{feat['id']}-{i}"
    else:
        # coding features all have same ID but are enumerated
        for sep in ['.', '-', '_']:
            if len(set([sep.join(id.split(sep)[:-1]) for id in cf_ids])) == 1:
                unique_coding_id = sep.join(cf_ids[0].split(sep)[:-1])
                gene['unique_coding_id'] = unique_coding_id
                break

    if not longest_transcript:
        if unique_coding_id:
            # put unambigous id of coding features as transcript id
            gene['transcript_id'] = unique_coding_id
        else:
            # put gene id as transcript id in case coding features dont have unambiguous id
            gene['transcript_id'] = gene.get('id')
    else:
        gene['transcript_id'] = longest_transcript.get('id')

    gene['coding_features'] = [feat.get('id') for feat in coding_features]
    gene['coding_coordinates'] = tuple([(feat.get('start'),feat.get('end')) for feat in coding_features])
    gene['gene_biotype'] = gene_biotype

    for feat in coding_features:
        feat['gene_id'] = gene.get('id')

    return [gene] + coding_features


def parse_genes(cfg):
    """ Extract genes and their transcript type from the GFF """

    feature_dicts = []
    gene = None

    with open(cfg.gff_path, 'r') as gff_file:
        for line_num, line in enumerate(gff_file):
            #print(line_num)

            # skip comments
            if line.startswith('#'):
                # stop parsing when FASTA block is reached
                if line.startswith('##FASTA'):
                    break
                else:
                    continue


            spline = line.strip().split('\t')
            if len(spline) == 9:
                feature_dict = spline2dict(spline, cfg.include_pseudogenes)
            else:
                spline.append(f"\tID={spline[2]}-{line_num}")
                feature_dict = spline2dict(spline, cfg.include_pseudogenes)

            if feature_dict.get('type') == 'gene':
                if gene:
                    feature_dicts += process_gene(gene, parent_ids, transcripts, gene_features, non_transcript_types)

                gene = feature_dict
                parent_ids = [gene.get('id')]
                transcripts = [] # list of all features with gene=Parent
                gene_features = [] # list of all other features connected to gene
                non_transcript_types = set()
            elif gene and feature_dict.get('parent_id') == gene.get('id'):
                transcripts.append(feature_dict)
                parent_ids.append(feature_dict.get('id'))
            elif gene and feature_dict.get('parent_id') in parent_ids:
                gene_features.append(feature_dict)
                parent_ids.append(feature_dict.get('id'))
                non_transcript_types.add(feature_dict.get('type'))
            else: # feature is not new gene but is also not connected to other gene features
                if gene:
                    feature_dicts += process_gene(gene, parent_ids, transcripts, gene_features, non_transcript_types)
                gene = None
        if gene:
            feature_dicts += process_gene(gene, parent_ids, transcripts,
                                          gene_features, non_transcript_types)

    gff_df = pd.DataFrame.from_dict(feature_dicts)

    return gff_df


def set_gene_neighbours(gff_df):


    # for scaffold_id in gff_df['scaffold'].unique():
    #     genes = gff_df.loc[(gff_df['scaffold'] == scaffold_id) & (gff_df['type'] == 'gene')].sort_values(['start'])
    #     gene_order = genes.index.to_list()
    #     gff_df.loc[gene_order,'upstream_gene'] = ([None] + gene_order)[:-1]
    #     gff_df.loc[gene_order,'downstream_gene'] = (gene_order + [None])[1:]

    cur_gene, cur_contig = None, None
    upstream_gene = None
    for row in gff_df.loc[gff_df['type'] == 'gene'].itertuples():
        #print(row)
        if cur_contig == row.scaffold:
            gff_df.at[cur_gene, 'upstream_gene'] = upstream_gene
            gff_df.at[cur_gene, 'downstream_gene'] = row.Index
            upstream_gene = cur_gene
            cur_gene = row.Index
        else:
            if cur_gene:
                gff_df.at[cur_gene, 'upstream_gene'] = upstream_gene
                gff_df.at[cur_gene, 'downstream_gene'] = None

            upstream_gene = None
            cur_contig = row.scaffold
            cur_gene = row.Index

    gff_df.at[cur_gene, 'upstream_gene'] = upstream_gene
    gff_df.at[cur_gene, 'downstream_gene'] = None

    return gff_df


def match_fasta2id(cfg, gff_df):
    """ Match the header of a precomputed protein FASTA file to gene ID

    :param geneid_in_fasta:
    :return:
    """

    header_dict = {}
    if cfg.extract_proteins:
        gff_df['fasta_header'] = gff_df['transcript_id']
    else:
        if cfg.prot2gene_mapper:
            # file that maps protein fasta file headers to gene IDs
            #TODO: change the delimiter
            prot2gene = pd.read_csv(cfg.prot2gene_mapper, sep='\t',
                                    names=['fasta_header', 'g_name'],
                                    index_col='g_name')

            gff_df = gff_df.join(prot2gene)

            unmapped_header = list(set(prot2gene['fasta_header']) - set(gff_df['fasta_header']))

        else:
            # if proteins will not be extracted by taxaminer
            # check of which feature type the ID is used in the protein FASTA
            with open(cfg.proteins_path, 'r') as p_file:
                for line in p_file:
                    if line.startswith('>'):
                        protein_id = line.strip().split()[0].split('|')[0][1:]
                        # diamond uses everything until first whitespace as ID
                        header_dict[protein_id] = line[1:].strip().split()[0]
                    else:
                        continue
            # if protein was passed by user, 'transcript_id' resembles the id
            # of the feature that denotes the header in the fasta file
            gff_df['fasta_header'] = gff_df['transcript_id'].map(header_dict)
            gff_df.loc[(gff_df['type'] == 'gene') & (gff_df['fasta_header'].isna()),'fasta_header'] = gff_df['unique_coding_id'].map(header_dict)

            unmapped_header = [key for key in header_dict.keys() if
                               key not in gff_df['fasta_header']]

        # check for protein coding genes that have no assigned protein fasta header
        missing_header = gff_df['fasta_header'].isna()
        # non-gene features and non-protein-coding genes don't require fasta header
        missing_header.loc[gff_df['type'] != 'gene'] = False
        missing_header.loc[(gff_df['type'] == 'gene') & (~gff_df['gene_biotype'].isin(['mRNA', 'exon', 'CDS']))] = False
        if missing_header.any():
            if gff_df[missing_header].shape[0] == gff_df.loc[(gff_df['type'] == 'gene') & (gff_df['gene_biotype'].isin(['mRNA', 'exon', 'CDS']))].shape[0]:
                # not any gene was mapped to a fasta header
                logging.error(
                    f"Unable to map the gene IDs to the headers in the protein "
                    f"FASTA file. Header must start with ID of either the gene, "
                    f"the transcript or the coding sequence in the GFF. "
                    f"If mapping continues to fail, make use of mapping file. "
                    f"For details see documentation")
                sys.exit(1)
            else:
                if not unmapped_header:
                    # if every fasta file header was mapped
                    logging.info(f"All headers of the fasta file were mapped to "
                                 f"a gene ID. Additional {gff_df.loc[gff_df['type'] == 'gene'][missing_header].shape[0]} "
                                 f"genes remain without a protein sequence.")
                else:
                    logging.info(
                        f"Unable to map the following {gff_df[missing_header].shape[0]} "
                        f"gene IDs to a header in the protein FASTA file:\n"
                        f"{gff_df[missing_header].index.to_list()}")
        else:
            logging.info(f"All genes attributed with protein sequence header.")


    gff_df['diamond_header'] = [str(x).split()[0] for x in gff_df['fasta_header'].tolist()]

    return gff_df


def process(cfg):

    #logging.debug(">> checking gff")
    # TODO: rewrite GFF check
    logging.debug(">> parsing file")
    gff_df = parse_genes(cfg)
    gff_df = strip_ids(gff_df)
    gff_df = add_missing_ids(gff_df)
    gff_df = remove_duplicate_ids(gff_df)
    gff_df = set_ids2index(gff_df)
    logging.debug(">> assessing gene neighbours")
    gff_df = set_gene_neighbours(gff_df)
    logging.debug(">> matching protein fasta headers to gene ids")
    gff_df = match_fasta2id(cfg, gff_df)

    return gff_df


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    process(cfg)


if __name__ == '__main__':
    main()
