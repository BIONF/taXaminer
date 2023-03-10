#!/usr/bin/env python

"""Prepare data for running taXaminer

#TODO: parse gff file here
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
from Bio import Seq
import pandas as pd

from . import checkInput

class Feature:
    """Object for GFF feature
    """
    def __init__(self, contig, info_dict):

        self.contig = contig
        self.id = info_dict.get('ID')
        self.parent = info_dict.get('Parent')
        self.start = info_dict.get('start')
        self.end = info_dict.get('end')
        self.strand = info_dict.get('strand')
        self.phase = info_dict.get('phase')
        self.biotype = info_dict.get('biotype') if info_dict.get('biotype') else info_dict.get('gene_biotype')
        self.transl_table = info_dict.get('transl_table')

        self.children = {}
        self.parsed_parent = None

        # gene only
        self.transcripts = {} # key = transcript ID, value = CDS length
        self.cdss = {} # key = transcript ID, value = [CDS IDs]
        self.coordinates = [] # coordinates of CDS of currently longest transcript
        self.phased_coordinates = []    # coordinates of CDS of currently longest transcript
                                        # adjusted for the information of phase attribute
        self.prot_seqs = {'phased':'', 'non-phased': ''} # store both seqs



def get_attr(attr_list, attr, na_return=None):
    """
    :param attr_list: raw GFF attribute list
    :param attr: attribute to return
    :param na_return: value to return if attribute does not exists
    :return: value of attribute if exists, else na_return
    """

    if attr in attr_list:
        return attr_list.split(f'{attr}=')[1].split(';')[0]
    else:
        return na_return


def parse_file(cfg):

    gff_df = pd.read_csv(cfg.gff_path, sep='\t', names=['scaffold', 'source', 'type',
                                               'start', 'end', 'score', 'strand', 'add_attrs'])

    gff_df['id'] = gff_df.apply(
        lambda row: get_attr(row['add_attrs'],'ID'), axis=1)
    gff_df['parent'] = gff_df.apply(
        lambda row: get_attr(row['add_attrs'],'Parent'), axis=1)
    gff_df['transl_table'] = gff_df.apply(
        lambda row: get_attr(row['add_attrs'], 'transl_table', 1), axis=1)


    return gff_df


def process(cfg):

    gff_df = parse_file(cfg)



def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    process(cfg)


if __name__ == '__main__':
    main()
