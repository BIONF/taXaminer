#!/usr/bin/env python

"""Perform taxonomic assignment

Expects processed config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

from . import checkInput
from . import compTaxonomicAssignment

import taxopy
import sys
import csv
import pathlib
import subprocess
import logging
import pandas as pd

def percentage_target(contig_genes):

    target_count = contig_genes.value_counts("is_target")
    return target_count[1]/target_count.values.sum()


def contig_is_target(contig_genes):

    if not 0 in contig_genes.is_target.values:
        return 1
    else:
        return 0


def comp_contig_lca(contig_genes):
    """


    :return:
    """

    gene_assignments = contig_genes['taxon_assignmentID'].unique()
    lca = compTaxonomicAssignment.compute_lca(gene_assignments)
    return lca


###############################################################################
###############################################################################


def process_assignments(cfg, gff_df, all_data_df, TAX_DB_local):
    """

    Args:
      cfg:

    Returns:

    """

    global TAX_DB
    TAX_DB = TAX_DB_local

    contigs = pd.DataFrame(index=all_data_df.c_name.unique())

    for contig in contigs:
        contig_genes = all_data_df.loc[all_data_df['c_name'] == contig, :]

        contigs.at[contig, 'assignments_lca'] = comp_contig_lca(contig_genes)
        contigs.at[contig, 'is_target'] = contig_is_target(contig_genes)
        contigs.at[contig, 'percentage_target'] = percentage_target(contig_genes)
        contigs.at[contig, 'lca'] = comp_contig_lca(contig_genes)

    return



def main():
    """ """
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    process_assignments(cfg, gff_df, taxonomic_assignment, TAX_DB_local)


if __name__ == '__main__':
    main()
