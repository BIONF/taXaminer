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
__version__ = "0.6.0"

import logging
import sys
import pathlib
import subprocess
import time
import timeit

from Bio import Seq
import pandas as pd

from . import checkInput

def check_gff(cfg):

    gene, transcript, exon, cds = None, None, None, None
    features = {}
    with open(cfg.gff_path, 'r') as gff_file:
        for line in gff_file:
            # skip comments
            if line.startswith('#'):
                continue
            spline = line.strip().split('\t')
            feature_dict = spline2dict(spline)
            if feature_dict.get('type') == 'gene':
                # read only first gene
                if gene:
                    break
                gene = feature_dict
                features[gene.get('id')] = gene
            elif feature_dict.get('type') == 'mRNA':
                transcript = feature_dict
                features[transcript.get('id')] = transcript
            elif feature_dict.get('type') == 'exon':
                exon = feature_dict
                features[exon.get('id')] = exon
            elif feature_dict.get('type') == 'CDS':
                cds = feature_dict
                features[cds.get('id')] = cds
            else:
                continue

    transcript_type = 'gene' # feature type of which to use in the protein FASTA
    if transcript:
        transcript_avail = True
        transcript_type = 'mRNA'
    if exon:
        exon_avail = True
        if transcript_type == 'gene':
            transcript_type = 'exon'
    if cds:
        cds_avail = True
        coding_feature = cds
        if transcript_type == 'gene' or transcript_type == 'exon':
            transcript_type = 'CDS'
    else:
        if exon:
            coding_feature = exon
        else:
            if transcript:
                coding_feature = transcript
            else:
                coding_feature = gene

    coding_type = coding_feature.get('type')
    parent_path = [coding_type]
    parent = coding_feature.get('id')
    while not features.get(parent).get('type') == 'gene':
        parent = features.get(parent).get('parent_id')
        parent_path.append(features.get(parent).get('type'))

    transcript_type = None
    if not cfg.extract_proteins:
        # if proteins will not be extracted by taxaminer
        # check of which feature type the ID is used in the protein FASTA
        with open(cfg.proteins_path, 'r') as p_file:
            for line in p_file:
                if line.startswith('>'):
                    protein_id = line.strip().split('|')[0][1:]
                    for id, feat in features.items():
                        if id == protein_id:
                            transcript_type = feat.get('type')
                            break
                        elif feat.get('add_attrs').get('ID') == protein_id:
                            transcript_type = feat.get('type')
                            break
                else:
                    if transcript_type:
                        break
                    continue

    return coding_type, transcript_type, parent_path


def match_fasta2id(cfg, gff_df):
    """ Match the header of a precomputed protein FASTA file to gene ID

    :return:
    """

    header_dict = {}
    if not cfg.extract_proteins:
        # if proteins will not be extracted by taxaminer
        # check of which feature type the ID is used in the protein FASTA
        with open(cfg.proteins_path, 'r') as p_file:
            for line in p_file:
                if line.startswith('>'):
                    protein_id = line.strip().split('|')[0][1:]
                    # diamond uses everything until first whitespace as ID
                    header_dict[protein_id] = line[1:].strip().split()[0]
                else:
                    continue
        gff_df['fasta_header'] = gff_df['transcript_id'].map(header_dict)
    else:
        gff_df['fasta_header'] = gff_df['transcript_id']


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
    if not id:
        return None

    for prefix in ['gene:', 'gene-', 'transcript:', 'transcript-',
                   'rna:', 'rna-', 'cds:', 'cds-']:
        id = remove_prefix(id, prefix)
    return id


def attrs2dict(attr_list):
    """ Parse the GFF additional attributes list into a dict

    :param attr_list:
    :return:
    """

    return {attr.split('=')[0]: attr.split('=')[1] for attr in
            attr_list.split(';')}


def spline2dict(spline):

    spline_dict = {}
    # direct inference
    spline_dict['scaffold'] = spline[0]
    spline_dict['source'] = spline[1]
    spline_dict['type'] = spline[2]
    spline_dict['start'] = int(spline[3])
    spline_dict['end'] = int(spline[4])
    spline_dict['score'] = spline[5]
    spline_dict['strand'] = spline[6]
    spline_dict['phase'] = spline[7]
    spline_dict['add_attrs'] = attrs2dict(spline[8])
    spline_dict['add_attrs']['ID'] = strip_ID(spline_dict['add_attrs'].get('ID'))
    spline_dict['add_attrs']['Parent'] = strip_ID(spline_dict['add_attrs'].get('Parent'))
    spline_dict['id'] = spline_dict['add_attrs'].get('ID')
    spline_dict['parent_id'] = strip_ID(spline_dict['add_attrs'].get('Parent'))
    transl_table = spline_dict['add_attrs'].get('transl_table')
    spline_dict['transl_table'] = transl_table if transl_table else '1'

    # indirect inference
    spline_dict['gene_id'] = spline_dict['id'] if spline_dict['type'] == 'gene' else None
    spline_dict['length'] = spline_dict['end'] - spline_dict['start'] + 1
    if not spline_dict.get('id'):
        # no ID information only acceptable for features without subfeatures
        spline_dict['id'] = f"{spline_dict.get('parent_id')}-{spline_dict['type']}-{spline_dict['start']}"

    return spline_dict


def replace_longest_transcript(gene, transcript, transcript_type, transcript_cds_features):

    cds_ids = [cds.get('id') for cds in transcript_cds_features]
    if len(set(cds_ids)) != len(cds_ids):
        cds_ids = []
        for i, cds in enumerate(transcript_cds_features):
            cds['id'] = f"{cds.get('id')}-{i}"
            cds_ids.append(cds.get('id'))
    cds_coordindates = [(cds.get('start'), cds.get('end')) for cds in
                        transcript_cds_features]
    transcript['coding_features'] = cds_ids
    gene['coding_features'] = cds_ids
    transcript['coding_coordindates'] = cds_coordindates
    gene['coding_coordindates'] = cds_coordindates
    if transcript_type:
        if transcript_type == transcript.get('type'):
            gene['transcript_id'] = transcript.get('id')
        else:
            gene['transcript_id'] = transcript_cds_features[0].get('add_attrs').get('ID')
    else:
        gene['transcript_id'] = transcript.get('id')
    # use translational table of CDS feature is gene does not have info
    if gene.get('transl_table') == '1' and transcript_cds_features[0].get('transl_table') != '1':
        # if coding feature has different transl table info than default use this
        gene['transl_table'] = transcript_cds_features[0].get('transl_table')
    gene_cds_features = [transcript]
    gene_cds_features += transcript_cds_features

    return gene_cds_features


def save_gene_features(feature_dict, pandas_row_list, transcript, gene,
                       transcript_cds_features, transcript_cds_length,
                       gene_cds_features, max_cds_len, upstream_gene, transcript_type):
    # process previous gene
    if transcript:
        if transcript_cds_length > max_cds_len:
            gene_cds_features = replace_longest_transcript(
                gene, transcript, transcript_type, transcript_cds_features)
    if gene:
        # add gene neighbours
        downstream_gene = feature_dict.get('id')
        gene['upstream_gene'] = upstream_gene
        gene['downstream_gene'] = downstream_gene
        upstream_gene = gene.get('id')
        if gene.get('coding_features'):
            pandas_row_list.append(gene)
            pandas_row_list += gene_cds_features
        else:
            # gene has no subfeatures -> gene = transcript = CDS
            gene['transcript_id'] = gene.get('id')
            gene['coding_features'] = [gene.get('id')]
            pandas_row_list.append(gene)

    return upstream_gene


def parse_file(cfg, coding_type, transcript_type, parent_path):
    """

    :param coding_type:
    :param parent_path:
    :param cfg:
    :return:
    """

    pandas_row_list = []

    gene, transcript, exon, cds = None, None, None, None
    upstream_gene = None
    transcript_cds_features, transcript_cds_length, gene_cds_features = None, None, None
    max_cds_len = None

    with open(cfg.gff_path, 'r') as gff_file:
        for line in gff_file:

            # skip comments
            if line.startswith('#'):
                # stop parsing when FASTA block is reached
                if line.startswith('##FASTA'):
                    break
                else:
                    continue

            spline = line.strip().split('\t')
            if not spline[2] in parent_path:
                continue
            feature_dict = spline2dict(spline)

            # feature == gene
            if feature_dict.get('type') == parent_path[-1]:
                # process previous gene
                upstream_gene = save_gene_features(
                    feature_dict, pandas_row_list, transcript, gene,
                    transcript_cds_features, transcript_cds_length,
                    gene_cds_features, max_cds_len, upstream_gene,
                    transcript_type)

                # init new gene
                gene = feature_dict
                max_cds_len = 0
                gene_cds_features = []
                transcript, cds = None, None

            # feature == coding feature
            elif feature_dict.get('type') == parent_path[0]:
                cds = feature_dict
                cds['gene_id'] = gene.get('id')
                if cds.get('gene_id') == cds.get('parent_id'):
                    if cds.get('length') > max_cds_len:
                        max_cds_len = cds.get('length')
                        gene_cds_features = [cds]
                        # update gene with info from cds
                        if transcript_type:
                            if transcript_type == gene.get('type'):
                                gene['transcript_id'] = gene.get('id')
                            else:
                                gene['transcript_id'] = cds.get('id')
                        else:
                            gene['transcript_id'] = gene.get('id')
                        gene['coding_features'] = [cds.get('id')]
                        gene['coding_coordindates'] = [(cds.get('start'), cds.get('end'))]
                        # use translational table of CDS feature is gene has default info
                        if gene.get('transl_table') == '1' and cds.get('transl_table') != '1':
                            gene['transl_table'] = cds.get('transl_table')

                else:
                    # add to list of CDS for transcript
                    transcript_cds_features.append(cds)
                    transcript_cds_length += cds.get('length')

             # feature between gene and CDS
            else:
                if feature_dict.get('type') == 'mRNA':
                    # process previous transcript
                    if transcript:
                        if transcript_cds_length > max_cds_len:
                            max_cds_len = transcript_cds_length
                            gene_cds_features = replace_longest_transcript(
                                gene, transcript, transcript_type,
                                transcript_cds_features)
                    # init new transcript
                    transcript = feature_dict
                    transcript['gene_id'] = gene.get('id')
                    transcript_cds_features = []
                    transcript_cds_length = 0
                else:
                    # feature only required for matching to gene ID, which
                    # is already fulfilled by order of the GFF file
                    continue

    save_gene_features(
        {}, pandas_row_list, transcript, gene, transcript_cds_features,
        transcript_cds_length, gene_cds_features, max_cds_len, upstream_gene,
        transcript_type)

    gff_df = pd.DataFrame(pandas_row_list)
    # if features don't have ID attribute, concatenate type and line number
    no_id = gff_df['id'].isnull()
    gff_df.loc[no_id, 'id'] = gff_df['type'] + '-' + gff_df.index.astype(str)
    # append start coordinate to duplicate ids
    duplicate_index = gff_df.duplicated(subset='id', keep=False)
    gff_df.loc[duplicate_index, 'id'] = gff_df['id'] + '-' + gff_df[
        'start'].astype(str)

    gff_df.set_index('id', inplace=True)

    return gff_df


def process(cfg):
    coding_type, transcript_type, parent_path = check_gff(cfg)
    gff_df = parse_file(cfg, coding_type, transcript_type, parent_path)
    match_fasta2id(cfg, gff_df)

    return gff_df


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    process(cfg)


if __name__ == '__main__':
    main()
