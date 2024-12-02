#!/usr/bin/env python

"""Perform taxonomic assignment

Expects processed config file
"""
__author__ = "Freya Arthen"

import time

from tqdm import tqdm

from . import checkInput
from . import compTaxonomicAssignment

import sys
import pandas as pd
import multiprocessing as mp

def chunks(lst, n):
    """Yield n successive chunks from lst."""
    size = (len(lst) // n) + 1
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

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


def comp_contig_lca(contig_genes, missing_taxids, TAX_DB):
    """

    :param TAX_DB:
    :param missing_ids:
    :return:
    """

    gene_assignments = contig_genes['taxon_assignmentID'].dropna().unique()
    if len(gene_assignments) == 0:
        return 'Unassigned'
    lca = compTaxonomicAssignment.compute_lca(gene_assignments, missing_taxids, TAX_DB)
    if lca:
        return lca
    else:
        return None


def comp_majority_assignment(contig_genes, fraction, missing_taxids, TAX_DB):


    gene_assignments = contig_genes['taxon_assignmentID'].dropna().unique()
    if len(gene_assignments) == 0:
        return 'Unassigned'
    majority = compTaxonomicAssignment.compute_majority_taxon(gene_assignments,
                                                              fraction, missing_taxids, TAX_DB)
    return majority.name


def comp_most_abundant_taxon(contig_genes):

    candidate_taxon = contig_genes['taxon_assignment'].dropna().value_counts().index
    if candidate_taxon.empty:
        return None
    else:
        return candidate_taxon[0]


def monitor_coverage(contig_genes):

    covs = contig_genes.filter(regex=("(g_cov_[0-9]*)|is_target"))
    target_average = covs.loc[covs['is_target'] == 1].drop(columns='is_target').mean()
    target_std = covs.loc[covs['is_target'] == 1].drop(columns='is_target').std()
    cont_average = covs.loc[covs['is_target'] == 0].drop(columns='is_target').mean()
    cont_std = covs.loc[covs['is_target'] == 0].drop(columns='is_target').std()


def compute_contig_stats(params):

    (chunk, all_data_df, missing_taxids, TAX_DB) = params

    df_list = []

    contigs_df_header = ['c_name', 'num_of_genes', 'percentage_target', 'lca',
        'most_abundant_taxon', 'target_bitscore_mean', 'target_bitscore_sd',
        'target_evalue_mean', 'target_evalue_sd', 'other_bitscore_mean',
        'other_bitscore_sd', 'other_evalue_mean',  'other_evalue_sd']

    for contig_id in chunk:

        contig_genes = all_data_df.loc[all_data_df['c_name'] == contig_id, :]
        contig_row = [contig_id, contig_genes.iloc[0].c_num_of_genes]

        if contig_genes.empty:
            df_list.append(contig_row)
            continue
        contig_lca = comp_contig_lca(contig_genes, missing_taxids, TAX_DB)
        if not contig_lca or contig_lca == "Unassigned":
            contig_row.append("Unassigned")
            contig_row.append(None)

        else:
            contig_row.append(contig_lca.name)
            contig_row.append(contig_lca.taxid)

        contig_row.append(percentage_target(contig_genes))
        contig_row.append(comp_most_abundant_taxon(contig_genes))

        target_genes = contig_genes.loc[contig_genes['is_target'] == 1]
        non_target_genes = contig_genes.loc[contig_genes['is_target'] == 0]
        # target_bitscore_mean
        contig_row.append(round(target_genes.bh_bitscore.astype(float).mean(), 2))
        # target_bitscore_sd
        contig_row.append(round(target_genes.bh_bitscore.astype(float).std(), 2))
        # target_evalue_mean
        contig_row.append(target_genes.bh_evalue.astype(float).mean())
        # target_evalue_sd
        contig_row.append(target_genes.bh_evalue.astype(float).std())
        # other_bitscore_mean
        contig_row.append(round(non_target_genes.bh_bitscore.astype(float).mean(), 2))
        # other_bitscore_sd
        contig_row.append(round(non_target_genes.bh_bitscore.astype(float).std(), 2))
        # other_evalue_mean
        contig_row.append(non_target_genes.bh_evalue.astype(float).mean())
        # other_evalue_sd
        contig_row.append(non_target_genes.bh_evalue.astype(float).std())

        df_list.append(contig_row)


    return df_list, missing_taxids

###############################################################################
###############################################################################


def process_assignments(cfg, gff_df, all_data_df, TAX_DB):
    """

    Args:
      cfg:

    Returns:

    """

    pool = mp.Pool(cfg.threads)

    missing_taxids = set()

    contigs_df_header = ['c_name', 'num_of_genes', 'lca', 'lcaID', 'percentage_target',
        'most_abundant_taxon', 'target_bitscore_mean', 'target_bitscore_sd',
        'target_evalue_mean', 'target_evalue_sd', 'other_bitscore_mean',
        'other_bitscore_sd', 'other_evalue_mean',  'other_evalue_sd']

    contig_ids = all_data_df.c_name.unique()

    chunked_contig_ids = list(chunks(contig_ids, cfg.threads))

    df_row_list = [contigs_df_header]
    for i in tqdm(pool.imap_unordered(compute_contig_stats,
            [(chunk, all_data_df, missing_taxids, TAX_DB) for
             chunk in chunked_contig_ids]),
            total=len(chunked_contig_ids),
            bar_format='{l_bar}{bar}| [{elapsed}<{remaining}, ' '{rate_fmt}{postfix}]'):
        # add data to df list
        df_row_list += i[0]
        missing_taxids = missing_taxids | i[1]

    contigs = pd.DataFrame(df_row_list)
    contigs.to_csv(f"{cfg.output_path}taxonomic_assignment/contig_assignments.csv", index=False, header=False)


def main():
    """ """
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    process_assignments(cfg, gff_df, taxonomic_assignment, TAX_DB)


if __name__ == '__main__':
    main()
