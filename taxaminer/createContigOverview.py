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
    if 1 in target_count:
        return target_count[1]/target_count.values.sum()
    else:
        return 0


def contig_is_target(contig_genes):

    if not 0 in contig_genes.is_target.values:
        return 1
    else:
        return 0


def comp_contig_lca(contig_genes):
    """


    :return:
    """

    gene_assignments = contig_genes['taxon_assignmentID'].dropna().unique()
    lca = compTaxonomicAssignment.compute_lca(gene_assignments)
    return lca.name

def comp_majority_assignment(contig_genes, fraction):


    gene_assignments = contig_genes['taxon_assignmentID'].dropna().unique()
    majority = compTaxonomicAssignment.compute_majority_taxon(
        gene_assignments, fraction)
    return majority.name

def monitor_coverage(contig_genes):

    covs = contig_genes.filter(regex=("(g_cov_[0-9]*)|is_target"))
    target_average = covs.loc[covs['is_target'] == 1].drop(columns='is_target').mean()
    target_std = covs.loc[covs['is_target'] == 1].drop(columns='is_target').std()
    cont_average = covs.loc[covs['is_target'] == 0].drop(columns='is_target').mean()
    cont_std = covs.loc[covs['is_target'] == 0].drop(columns='is_target').std()


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

    contig_ids = all_data_df.c_name.unique()
    contig_dict = {
        'num_of_genes': [None for i in range(contig_ids.size)],
        'percentage_target': [None for i in range(contig_ids.size)],
        'lca': [None for i in range(contig_ids.size)],
        'most_abundant_taxon': [None for i in range(contig_ids.size)]
    }
    contigs = pd.DataFrame(contig_dict, index=contig_ids)
    contigs['num_of_genes'] = contigs.index.map(
        all_data_df[['c_name', 'c_num_of_genes']].set_index('c_name').to_dict().get('c_num_of_genes'))

    for contig in contigs.itertuples():
        contig_genes = all_data_df.loc[all_data_df['c_name'] == contig.Index, :]

        if contig_genes.empty:
            continue
        contigs.at[contig.Index, 'lca'] = comp_contig_lca(contig_genes)
        contigs.at[contig.Index, 'percentage_target'] = percentage_target(
            contig_genes)
        #contigs.at[contig.Index, 'is_target'] = contig_is_target(contig_genes)
        # add taxon which represents x% of the data
        #contigs.at[contig.Index, 'majority_75_taxon'] = comp_majority_assignment(contig_genes, 0.75)
        contigs.at[
            contig.Index, 'most_abundant_taxon'] = contig_genes['taxon_assignment'].dropna().value_counts().index[0]
        #contigs.at[contig.Index, 'all_assignments'] = ';'.join(contig_genes['taxon_assignment'].dropna().unique().tolist())

        # if cfg.include_coverage:
        #     contigs.at[
        #         contig.Index, 'putative_hgt'] = monitor_coverage(
        #         contig_genes)

    contigs.to_csv(f"{cfg.output_path}taxonomic_assignment/contig_assignments.csv", index_label='c_name')


def main():
    """ """
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    process_assignments(cfg, gff_df, taxonomic_assignment, TAX_DB_local)


if __name__ == '__main__':
    main()
