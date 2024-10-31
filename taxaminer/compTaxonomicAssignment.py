#!/usr/bin/env python

"""Perform taxonomic assignment

Expects processed config file
"""
__author__ = "Freya Arthen"


from . import checkInput

import taxopy
import sys
import pathlib
import subprocess
import logging
import pandas as pd
import multiprocessing as mp
from tqdm import tqdm
import psutil

###############################################################################
############################ HELPER FUNCTIONS #################################
###############################################################################

def get_id_for_rank_of_species(target_id, rank, missing_taxids, TAX_DB):
    """
    Returns ID of taxon at given rank for given taxon ID.

    Args:
      target_id(int): NCBI ID of taxon for wich to find ID at given rank
      rank(str): rank of which to get the NCBI ID of

    Returns:
       (int) NCBI ID of taxon at the given rank for given taxon ID
       :param TAX_DB:
     """
    query_taxon = init_taxopyon(target_id, missing_taxids, TAX_DB)
    query_ranks = query_taxon.rank_name_dictionary
    rank_of_interest = query_ranks.get(rank)
    # get ID for the specific class (function taxopy.taxid_from_name() might get multiple IDs because of homonyms)
    if rank_of_interest:
        rank_index = query_taxon.name_lineage.index(rank_of_interest)
        rank_id = query_taxon.taxid_lineage[rank_index]
        return rank_id
    else:
        return None


def get_children_ids(TAX_DB, parent):
    """
    Returns list of taxon IDs which are children of given taxon

    Args:
      parent (int): NCBI taxon ID to get IDs of children for

    Returns:
        (list) list of taxon IDs which are children of given taxon
        :param TAX_DB:
    """
    # taxid2parent: dict with key = child, value = parent
    return [key for key in TAX_DB.taxid2parent.keys() if
            TAX_DB.taxid2parent[key] == parent]


def init_taxopyon(taxid, missing_taxids, TAX_DB):
    """ init taxopy taxon object
    :param TAX_DB:
    :param missing_taxids:
    """

    if ';LCA:' in str(taxid):
        taxid = taxid.split(';')[0]

    try:
        taxon = taxopy.Taxon(int(taxid), TAX_DB)
        return taxon
    except:
        try:
            new_taxid = TAX_DB.oldtaxid2newtaxid[int(taxid)]
            taxon = taxopy.Taxon(new_taxid, TAX_DB)
            return taxon
        except:
            missing_taxids.add(int(taxid))
            return None

###############################################################################
######################### PLOT LABEL COMPUTATION ##############################
###############################################################################

def get_children_for_label(taxa, labelID, missing_taxids, TAX_DB):
    """
    gets all descendants of label which are in lineage of one of the taxa the label represents
    also counts how many of these descendants are already in list of current labels

    Args:
      taxa:
      labelID:

    Returns:
    :param TAX_DB:

    """
    children = {}
    for taxID in taxa:
        taxon = init_taxopyon(taxID, missing_taxids, TAX_DB)
        child = taxon.taxid_lineage[max(taxon.taxid_lineage.index(labelID) - 1,
                                        0)]  # max 0 if taxon is end of lineage
        if child in children:
            children[child].add(taxID)
        else:
            children[child] = set([taxID])

    return children


def merge_assignments_at_rank(assignments_df, rank, missing_taxids, TAX_DB):
    """
    Merge taxonomic assignments at given rank.
    Args:
      genes(dict): dictionary of genes {g_name: gene object}
      rank(str): rank at which to merge all taxonomic assignments at
      :param TAX_DB:
    """

    for assignmentID in assignments_df['taxon_assignmentID'].dropna().unique():
        genes_w_assignment = (
                    assignments_df['taxon_assignmentID'] == assignmentID)
        taxon = init_taxopyon(assignmentID, missing_taxids, TAX_DB)
        ranks = taxon.rank_name_dictionary
        merger_rank = ranks.get(rank)
        if merger_rank:
            assignments_df.loc[genes_w_assignment, 'plot_label'] = merger_rank
            assignments_df.loc[genes_w_assignment, 'plot_labelID'] = taxon.rank_taxid_dictionary.get(rank)
        else:
            assignments_df.loc[genes_w_assignment, 'plot_label'] = assignments_df.loc[genes_w_assignment, 'taxon_assignment']
            assignments_df.loc[genes_w_assignment, 'plot_labelID'] = assignments_df.loc[genes_w_assignment, 'taxon_assignmentID']


def merge_assignments_at_id(assignments_df, ids, missing_taxids, TAX_DB):
    """
    Merge taxonomic assignments at given NCBI taxon IDs.

    Args:
      assignments_df(dict): dictionary of genes {g_name: gene object}
      ids(int,list): list or int of taxon ID(s) at which to merge
      :param TAX_DB:
    """

    if type(ids) != list:
        ids = [ids]

    for assignmentID in assignments_df['taxon_assignmentID'].dropna().unique():
        id_merger_taxon = None
        taxon = init_taxopyon(assignmentID, missing_taxids, TAX_DB)
        lineage = taxon.taxid_lineage
        for id in ids:
            if id in lineage:
                id_merger_taxon = init_taxopyon(id, missing_taxids, TAX_DB)
                lineage = id_merger_taxon.taxid_lineage
                break
        if id_merger_taxon:        
            genes_w_assignment = (assignments_df['taxon_assignmentID'] == assignmentID)
            assignments_df.loc[genes_w_assignment, 'plot_label'] = id_merger_taxon.name
            assignments_df.loc[genes_w_assignment, 'plot_labelID'] = id_merger_taxon.taxid

def rank_sorted_taxonids(ids, missing_taxids, TAX_DB):

    # use heuristic:
    # sort by length of the lineage (shortest first)
    # bacteria would be expanded first, since their lineages are often shorter
    sorted_ids = []
    delete_ids = []
    for rank in ['species', 'genus', 'family', 'order', 'class',
                 'phylum', 'kingdom', 'superkingdom']:
        for id in ids:
            taxon = init_taxopyon(id, missing_taxids, TAX_DB)
            if rank in taxon.rank_name_dictionary:
                sorted_ids.append(id)
                delete_ids.append(id)
        ids = [id for id in ids if not id in delete_ids]
        delete_ids = []
    sorted_ids += ids

    return(sorted_ids[::-1])



def iterative_merging_top(assignments_df, num_groups_plot, missing_taxids, TAX_DB):
    """
    identify at which rank in query lineage to merge best

    Args:
      genes_df:
      num_groups_plot:

    Returns:
    :param TAX_DB:
    :param missing_taxids:

    """

    unlabeled_genes = assignments_df.loc[assignments_df['plot_label'].isnull()]
    # taxa that have this label; label : set(taxa)
    labelID_taxa = {1: set(unlabeled_genes['taxon_assignmentID'].dropna().unique().astype(int))}

    update_dict = {'nonempty-dummy'}
    while update_dict:
        update_dict = {}
        del_dict_keys = set()
        remaining_labels_num = num_groups_plot - len(labelID_taxa)
        rank_sorted_labelIDs = rank_sorted_taxonids(labelID_taxa.keys(), missing_taxids, TAX_DB)
        for labelID in rank_sorted_labelIDs:
            taxa = labelID_taxa.get(labelID)
            if len(taxa) >= 2:
                lca_candidates = []
                for taxon in taxa:
                    lca_candidates.append(
                        init_taxopyon(taxon, missing_taxids, TAX_DB))
                lca = taxopy.find_lca(lca_candidates, TAX_DB)
                if lca.taxid != labelID:  # if it is possible to find a more specific label than the current one with the lca
                    del_dict_keys.add(labelID)
                    update_dict[lca.taxid] = taxa
                else:
                    # get all children for label that are in list taxa
                    children = get_children_for_label(taxa, labelID,
                                                      missing_taxids, TAX_DB)
                    existing_labels = len(children.keys()) - len(
                        list(set(set(children.keys()) - labelID_taxa.keys())))

                    # check if number of new labels (= num of children - labels that are already in lablels)
                    # fits in remaining number of possible new labels
                    if len(children) - existing_labels <= remaining_labels_num:
                        update_dict.update(children)
                        del_dict_keys.add(labelID)
                        remaining_labels_num -= (
                                    len(children) - existing_labels)
            else:
                pass
        for key in del_dict_keys:
            del labelID_taxa[key]
        labelID_taxa.update(update_dict)

    taxid2label = {}  # {taxonomic assignment : the label it gets}
    for labelID, taxa in labelID_taxa.items():
        if len(taxa) == 1:
            # min(taxa) instead of taxa[0] since taxa is a set
            # and therefor not subscriptable
            taxid2label[min(taxa)] = min(taxa)
        else:
            for taxon in taxa:
                taxid2label[taxon] = labelID

    assignments_df['plot_label'] = assignments_df.loc[~assignments_df['taxon_assignmentID'].isnull()].apply(
        lambda row: TAX_DB.taxid2name.get(taxid2label.get(int(row.taxon_assignmentID))) if (row.taxon_assignmentID and not row.plot_label) else row.plot_label, axis=1)
    assignments_df['plot_labelID'] = assignments_df.loc[~assignments_df['taxon_assignmentID'].isnull()].apply(
        lambda row: taxid2label.get(
            int(row.taxon_assignmentID)) if not row.plot_labelID else row.plot_labelID, axis=1)


    return assignments_df


def merge_labels(assignments_df, target_id, merging_labels, missing_taxids, TAX_DB):
    """
    Merge plot labels based on user input.

    Args:
      assignments_df(dict): dictionary of genes {g_name: gene object}
      target_id(int): NCBI taxon ID for query species
      merging_labels(list,int): user input for merging labels
        1. taxonomic rank with suffix '-all'
            -> merge all taxonomic assignment at that rank
        2. taxonomic rank
            -> merge query lineage at that rank
        3. list of NCBI taxon IDs
            -> merge at each of the taxa provided
            :param missing_taxids:
            :param TAX_DB:
    """

    merger_ids = []
    all_merger_rank = None
    query_merger_id = None

    if not merging_labels:
        merging_labels = []
    elif type(merging_labels) != list:
        if ',' in merging_labels:
            merging_labels = merging_labels.split(',')
        else:
            merging_labels = [merging_labels]

    for merger_label in merging_labels:

        if not merger_label:
            continue

        if type(merger_label) == str:
            # option 1: string end with '-all' -> every taxonomic assignment is merged at that rank
            # option 2: string is a NCBI Taxonomy ID in string format
            # option 3: string does not end with '-all' -> only taxonomic assignments with same rank <merger_label> as the query species are merged

            if merger_label.endswith('-all'):
                all_merger_rank = merger_label[:-len('-all')]
                break
            elif merger_label.isdigit():
                merger_ids.append(
                    int(merger_label))  # ncbiID of <merger_label> to merge at            
            else:
                # get ncbiID for rank <merger_label> of query
                query_merger_id = get_id_for_rank_of_species(target_id,
                                                             merger_label, missing_taxids, TAX_DB)
                break
        else:
            if type(merger_label) == int or merger_label.isdigit():
                merger_ids.append(
                    int(merger_label))  # ncbiID of <merger_label> to merge at

    # if there is an -all merger-label; this overwrites everything
    # else; merge at query_merger_id and merger_ids
    if all_merger_rank:
        merge_assignments_at_rank(assignments_df, all_merger_rank, missing_taxids, TAX_DB)
    else:
        if query_merger_id:
            merge_assignments_at_id(assignments_df, query_merger_id, missing_taxids, TAX_DB)
        if merger_ids:
            merge_assignments_at_id(assignments_df, merger_ids, missing_taxids, TAX_DB)

    return assignments_df


def set_unassigned_labels(assignments_df):
    """
    Set missing plot labels to taxonomic assignment or 'Unassigned'.

    Args:
      genes:
    """

    unassigned_genes = assignments_df['plot_labelID'].isnull()
    assignments_df.loc[unassigned_genes, 'plot_label'] = 'Unassigned'
    return assignments_df


def ident_query_label(assignments_df, target_taxon, TAX_DB):
    """
    Identify which label represents the query species.

    Args:
      genes:
      target_id:

    Returns:
    :param TAX_DB:

    """

    target_lineage = target_taxon.taxid_lineage

    label_candidate = (None, len(target_lineage))
    for assignment in assignments_df['plot_labelID'].dropna().unique():
        if assignment in target_lineage:
            if target_lineage.index(assignment) < label_candidate[1]:
                label_candidate = (assignment, target_lineage.index(assignment))

    return TAX_DB.taxid2name[label_candidate[0]]


###############################################################################
########################## TAXONOMIC ASSIGNMENT ###############################
###############################################################################

def get_subset_id(target_id, subset_marker, missing_taxids, TAX_DB):
    """
    Return input if input.isdigit() or return ID for rank of query.

    Args:
      target_id:
      subset_marker:

    Returns:
    :param missing_taxids:
    :param TAX_DB:

    """
    if type(subset_marker) == int:
        return subset_marker
    elif subset_marker.isdigit():
        return int(subset_marker)
    else:
        return get_id_for_rank_of_species(target_id, subset_marker, missing_taxids, TAX_DB)


def get_total_and_current(text):
    """ Parse progress line of Diamond to get current status and total status

     Example texts: 'Processing query block 1', 'reference block 1/128', 'shape 4/16'

    :return:
    """
    both_values = text.split()[-1].split('/')
    if len(both_values) == 2:
        current_value, max_value = both_values
    else:
        current_value, max_value = both_values[0], both_values[0]
    return int(current_value), int(max_value)


def run_diamond(diamond_cmd):
    """
    perform DIAMOND run and catch errors

    Args:
      diamond_cmd:

    Returns:

    """

    logging.info('>> running DIAMOND')
    dmndOut = subprocess.Popen(diamond_cmd,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    init_pbar = False
    old_shapes = 0
    pbar = None
    full_log = []
    while dmndOut.poll() is None:
        log_line = dmndOut.stderr.readline().decode().strip()
        full_log.append(log_line)
        if 'reference block' in log_line:
            if not init_pbar:
                query_block, ref_block, shape = log_line[:-1].split(',')  # rm trailing .
                current_query_block, max_query_blocks = get_total_and_current(
                    query_block)
                current_ref_block, max_ref_blocks = get_total_and_current(
                    ref_block)
                current_shape, max_shapes = get_total_and_current(shape)
                total_shapes = max_query_blocks * max_ref_blocks * max_shapes
                pbar = tqdm(total=total_shapes,
                            bar_format='{l_bar}{bar}| [{elapsed}<{remaining}, ' '{rate_fmt}{postfix}]')
                init_pbar = True
            pbar.update(1)
        if log_line.startswith('Scoring parameters:'):
            logging.debug(log_line)
        elif log_line.startswith('Temporary directory:'):
            logging.debug(log_line)
        elif log_line.startswith(
                'Percentage range of top alignment score to report hits:'):
            logging.debug(log_line)
        elif log_line.startswith('Reference ='):
            logging.debug(log_line)
        elif log_line.startswith('Sequences ='):
            logging.debug(log_line)
        elif log_line.startswith('Letters ='):
            logging.debug(log_line)
        elif log_line.startswith('Block size ='):
            logging.debug(log_line)
        elif log_line.startswith('Total time ='):
            pbar.update(1)
            logging.info(f'DIAMOND runtime: {log_line.split(" = ")[1]}')
        elif ' pairwise alignments, ' in log_line:
            logging.info(log_line)
        elif ' queries aligned.' in log_line:
            logging.info(log_line)

    if dmndOut.returncode != 0:
        logging.error(f'Error running DIAMOND '
                      f'with following command:\n{" ".join(diamond_cmd)}')
        if any('Opening the database... No such file or directory' in line for line in full_log):
            logging.error('Database does not exists. Please check the taXaminer setup.')
        # TODO: check for output when putting empty protein file
        elif any('asnhdoiashd' in line for line in full_log):
            logging.error('Empty protein FASTA file. Please check your input.')
        else:
            log_string = "\n".join(full_log)
            logging.error(f'DIAMOND Error message:\n{log_string}')
        sys.exit()

    if pbar:
        pbar.close()


def perform_quick_search_1(cfg, perform_diamond, diamond_cmd,
                           quick_mode_search_rank, target_exclude, target_taxon, missing_taxids, TAX_DB):
    """
    Perform first DIAMOND run for quick assignment mode.

    Args:
      perform_diamond:
      diamond_cmd:
      tax_assignment_path_1:
      quick_mode_search_rank:
      target_exclude:
      target_id:
      :param missing_taxids:
      :param TAX_DB:
    """
    logging.info("Quick mode for taxonomic assignment selected")

    if perform_diamond:
        # quick_mode_search_rank can be either taxon ID or rank
        quick_mode_search_id = get_subset_id(target_taxon.taxid,
                                             quick_mode_search_rank,
                                             missing_taxids, TAX_DB)

        if target_exclude:
            q1_exclude = target_exclude + diamond_inclusion_by_exclusion(quick_mode_search_id, missing_taxids, TAX_DB)
        else:
            q1_exclude = ['--taxon-exclude'] + diamond_inclusion_by_exclusion(quick_mode_search_id,missing_taxids, TAX_DB)

        diamond_cmd_1 = diamond_cmd + ['-o', cfg.diamond_results_path] + q1_exclude
        logging.debug(f"Diamond command for "
                      f"quick search round 1:\n{diamond_cmd_1}")
        run_diamond(diamond_cmd_1)


def perform_quick_search_2(cfg, perform_diamond, diamond_cmd,
                           quick_mode_match_rank, assignment_df, tmp_prot_path,
                           target_exclude, target_taxon, missing_taxids, TAX_DB):
    """
    Perform second DIAMOND run for quick assignment mode.

    Args:
      perform_diamond:
      diamond_cmd:
      tax_assignment_path_2:
      quick_mode_match_rank:
      genes:
      tmp_prot_path:
      target_exclude:
      target_id:
      :param missing_taxids:
      :param TAX_DB:
    """

    if perform_diamond:
        quick_mode_match_id = get_subset_id(target_taxon.taxid,
                                            quick_mode_match_rank,
                                            missing_taxids, TAX_DB)

        filter_prots_quick_hits(assignment_df, quick_mode_match_id,
                                tmp_prot_path, tmp_prot_path, missing_taxids, TAX_DB)

        diamond_cmd_2 = diamond_cmd + ['-o', cfg.diamond_results_path] + target_exclude
        logging.debug(f"Diamond command for quick search round 2:\n{diamond_cmd_2}")
        run_diamond(diamond_cmd_2)


def perform_exhaustive_search(cfg, perform_diamond, diamond_cmd, target_exclude):
    """
    Perform  DIAMOND run for exhaustive assignment mode.

    Args:
      perform_diamond:
      diamond_cmd:
      tax_assignment_path:
      target_exclude:
    """
    logging.info("Exhaustive mode for taxonomic assignment selected")
    if perform_diamond:
        diamond_cmd += ['-o', cfg.diamond_results_path] + target_exclude
        logging.debug(
            f"Diamond command for exhaustive search:\n{diamond_cmd}")
        run_diamond(diamond_cmd)


def generalize_unclassified(assignments_df, missing_taxids, TAX_DB):

    # for assignment in tax assignemnts
    # if unclassified in taxonomy_string:
    #   move to last taxon without unclassified/environmental info

    unclassified_strings = ['unclassified', 'environmental', 'uncultured'] # incertae sedis?

    for assignmentID in assignments_df['taxon_assignmentID'].dropna().unique():
        new_assignment = None
        taxon = init_taxopyon(assignmentID, missing_taxids, TAX_DB)
        if any((string in taxon.name) for string in unclassified_strings):
            for taxon in zip(taxon.rank_name_dictionary.values(),
                             taxon.rank_taxid_dictionary.values()):
                if not any((string in taxon[0]) for string in unclassified_strings):
                    new_assignment = taxon
                    break
            if not new_assignment:
                print(taxon)
                continue
            genes_w_assignments = assignments_df['taxon_assignmentID'] == assignmentID
            assignments_df.loc[genes_w_assignments, 'taxon_assignment'] = new_assignment[0]
            assignments_df.loc[genes_w_assignments, 'taxon_assignmentID'] = new_assignment[1]


def assess_closest_of_hits(taxIDs, target_taxon, missing_taxids, TAX_DB):
    """
    for a list of taxa IDs, get the one where LCA with target_id is closest to query

    Args:
      taxIDs:
      target_id:

    Returns:
    :param TAX_DB:

    """

    target_lineage = target_taxon.taxid_lineage

    # (taxon, index in lineage, lca with target)
    closest_hit = (None, len(target_lineage) + 1)
    for taxID in taxIDs:     
        taxon = init_taxopyon(taxID, missing_taxids, TAX_DB)
        if taxon:
            lca = taxopy.find_lca([target_taxon, taxon], TAX_DB)
            lca_index = target_lineage.index(lca.taxid)
            # check if new best hit is closer to query in lineage as current best hit
            if closest_hit[1] > lca_index:              
                closest_hit = (taxon, lca_index, lca)

    return closest_hit[0], closest_hit[2]

def compute_majority_taxon(taxon_ids, fraction, missing_taxids, TAX_DB):
    """ compute the LCA of a list of taxon ID
    :param TAX_DB:
    """

    taxa_list = [init_taxopyon(taxon_id, missing_taxids, TAX_DB) for taxon_id in taxon_ids]
    nanfree_taxa_list = [taxon for taxon in taxa_list if taxon]
    if len(nanfree_taxa_list) > 1:
        return taxopy.find_majority_vote(nanfree_taxa_list, TAX_DB, fraction=fraction)
    else:
        return nanfree_taxa_list[0]


def compute_lca(taxon_ids, missing_taxids, TAX_DB):
    """ compute the LCA of a list of taxon ID
    :param TAX_DB:
    :param missing_taxids:
    """

    taxa_list = [init_taxopyon(taxon_id, missing_taxids, TAX_DB) for taxon_id in taxon_ids]
    nanfree_taxa_list = [taxon for taxon in taxa_list if taxon]

    if not nanfree_taxa_list:
        print(taxon_ids)
        return None
    if len(nanfree_taxa_list) > 1:
        return taxopy.find_lca(nanfree_taxa_list, TAX_DB)
    elif nanfree_taxa_list:
        return nanfree_taxa_list[0]
    else:
        return None


def write_hits2file(cfg, file_path, hits):

    with open(file_path, 'a') as file:
        for hit in hits:
            if cfg.slim_diamond_results:
                file.write('\t'.join(hit[:3]+hit[10:]) + '\n')
            else:
                file.write('\t'.join(hit) + '\n')


def process_multi_hit(hit, target_taxon, missing_taxids, TAX_DB):

    split_taxids = hit[-2].split(';')

    lca = compute_lca(split_taxids, missing_taxids, TAX_DB)
    if not lca:
        return ['', '']

    if lca.taxid in target_taxon.taxid_lineage:
        closest_hit, closest_lca = assess_closest_of_hits(split_taxids,
                                                          target_taxon,missing_taxids, TAX_DB)
        return [f"{closest_hit.taxid};LCA:{lca.taxid}",
                f"{closest_hit.name};LCA:{lca.name}"]
    else:
        return [f"{lca.taxid};LCA:{lca.taxid}", f"{lca.name};LCA:{lca.name}"]


def process_hit(hit, target_taxon, missing_taxids, TAX_DB):

    if ';' in hit[-2]:  # hit is a multihit
        cropped_hit = process_multi_hit(hit, target_taxon, missing_taxids,
                                        TAX_DB)
        return hit[:-2]+cropped_hit
    # elif 'N/A' in hit[-1]: # no taxon matched to accession number
    #     list.append(hit[:-2] + ['N/A', 'N/A'])
    else:
        return hit


def calc_assignment(params):
    """
    Calculate the LCA for each gene based on taxonomic hits.

    Args:
      genes:

    Returns:
    :param missing_taxids:
    :param TAX_DB:
    :param target_taxon:

    """
    hit_dict, assignments_df, target_taxon, missing_taxids, TAX_DB = params
    return_dict = {}
    unmapped_proteins = []


    for protein_name, raw_hitlist in hit_dict.items():

        gene = assignments_df.loc[assignments_df.diamond_header == protein_name]
        try:
            gene_name = gene.index.item()
        except:
            unmapped_proteins.append(protein_name)
            continue

        hitlist = [process_hit(hit.strip().split('\t'), target_taxon, missing_taxids, TAX_DB) for hit in raw_hitlist]
        hit_ids = [id[-2].split(';')[0] for id in hitlist if id[-2] != '']
        if not hit_ids:
            unmapped_proteins.append(protein_name)
            continue
        else:
            lca = compute_lca(hit_ids, missing_taxids, TAX_DB)
            if lca:
                if lca.taxid in target_taxon.taxid_lineage:
                    closest_hit, closest_lca =  assess_closest_of_hits(hit_ids,
                                                                      target_taxon, missing_taxids, TAX_DB)
                else:
                    closest_lca = None
            else:
                closest_lca = None
        # all hit IDs were unable to be recognized
        if not lca:
            lca = None

        return_dict[gene_name] = (hitlist, lca, closest_lca)

    return return_dict, unmapped_proteins


def add_ta2gene(gene_id, results, cfg, tax_assignment_path, assignments_df, target_taxon):

    (hitlist, lca, closest_lca) = results

    #gene_id = assignments_df.query(f'diamond_header == "{diamond_header}"').index.item()

    if not hitlist or not lca:
        return

    # lca is ubiquitous
    assignments_df.at[gene_id,'lcaID'] = lca.taxid
    assignments_df.at[gene_id,'lca'] = lca.name

    # best hit is ubiquitous
    assignments_df.at[gene_id, 'best_hitID'] = hitlist[0][-2].split(';')[0]
    assignments_df.at[gene_id, 'best_hit'] = hitlist[0][-1].split(';')[0]
    assignments_df.at[gene_id, 'bh_evalue'] = hitlist[0][10]
    assignments_df.at[gene_id, 'bh_pident'] = hitlist[0][2]
    assignments_df.at[gene_id, 'bh_bitscore'] = hitlist[0][11]

    # condition -> taxonomic assignment
    ## taxonomic assignment in target lineage -> refined LCA
    if lca.taxid in target_taxon.taxid_lineage:
        assignments_df.at[gene_id, 'refined_lcaID'] = closest_lca.taxid
        assignments_df.at[gene_id, 'refined_lca'] = closest_lca.name
        assignments_df.at[gene_id, 'taxon_assignmentID'] = closest_lca.taxid
        assignments_df.at[gene_id, 'taxon_assignment'] = closest_lca.name
        assignments_df.at[gene_id, 'is_target'] = 1

    ## taxonomic assignment in target genus -> target genus
    elif lca.rank_taxid_dictionary.get(cfg.ta_generalisation_rank) == target_taxon.rank_taxid_dictionary.get(cfg.ta_generalisation_rank):
        assignments_df.at[gene_id, 'taxon_assignmentID'] = lca.rank_taxid_dictionary.get(cfg.ta_generalisation_rank)
        assignments_df.at[gene_id, 'taxon_assignment'] = lca.rank_name_dictionary.get(cfg.ta_generalisation_rank)
        assignments_df.at[gene_id, 'is_target'] = 1

    ## none of the above -> LCA
    else:
        assignments_df.at[gene_id, 'taxon_assignmentID'] = lca.taxid
        assignments_df.at[gene_id, 'taxon_assignment'] = lca.name


    write_hits2file(cfg, tax_assignment_path, hitlist)


def read_hit_file(cfg, chunk_size):
    """
    Read DIAMOND output file and assign hits to genes.

    Args:
      tax_assignment_path:
      prots:
      target_id:

    Returns:

    """

    count = 0
    return_dict = {}
    with open(cfg.diamond_results_path, 'r') as diamond_hits:
        # tax assignments list is sorted by query sequence
        # thus, all hits for one gene are following one another

        first_line = next(diamond_hits)
        first_spline = first_line.strip().split('\t')
        if first_spline[0] == 'qseqid':
            #TODO: extend support to parsing already parsed diamond tables
            first_spline = next(diamond_hits).strip().split('\t')

        current_diamond_header = first_spline[0]
        gene_hits = [first_line]
        for line in diamond_hits:
            count += len(line)
            # first hit per gene is "best hit" in terms of alignment score
            if line.startswith(current_diamond_header): #diamond_header: #
                gene_hits.append(line)
            else:
                return_dict[current_diamond_header] = gene_hits
                if count >= chunk_size:
                    yield return_dict
                    return_dict = {}
                    count = 0

                # new gene
                current_diamond_header = line.split()[0]
                gene_hits = [line]

        return_dict[current_diamond_header] = gene_hits
        yield return_dict


def process_hits(cfg, tax_assignment_path, assignments_df, target_taxon, missing_taxids, TAX_DB):

    logging.info(f">> processing DIAMOND hits file")
    threads = max((cfg.threads // 4), 1)
    # start the processing
    pool = mp.Pool(threads)

    file_mem = pathlib.Path(cfg.diamond_results_path).stat().st_size
    available_mem = psutil.virtual_memory().available * 0.1
    chunk_size = available_mem / threads
    n_chunks = file_mem // chunk_size
    if (n_chunks < threads):
        n_chunks = threads
        chunk_size = file_mem // n_chunks

    # empty file
    with open(tax_assignment_path, 'w') as file:
        if cfg.slim_diamond_results:
            file.write('\t'.join(['qseqid', 'sseqid', 'pident', 'evalue',
                                  'bitscore', 'taxid', 'taxname']) + '\n')
        else:
            file.write('\t'.join(['qseqid', 'sseqid', 'pident', 'length',
                                  'mismatch', 'gapopen', 'qstart', 'qend',
                                  'sstart', 'send', 'evalue', 'bitscore',
                                  'taxid', 'taxname']) + '\n')

    chunks = read_hit_file(cfg, chunk_size)
    unmapped_proteins = []
    pbar = tqdm(total=n_chunks,
                bar_format='{l_bar}{bar}| [{elapsed}<{remaining}, ' '{rate_fmt}{postfix}]')

    for i in pool.imap_unordered(
            calc_assignment, [(chunk, assignments_df, target_taxon, missing_taxids, TAX_DB) for chunk in chunks]):
        pbar.update(1)
        for gene_name, results in i[0].items():
            add_ta2gene(gene_name, results, cfg, tax_assignment_path, assignments_df,
                        target_taxon)
        unmapped_proteins += i[1]


    if pbar.n != n_chunks:
        pbar.update(n_chunks-pbar.n)

    pbar.close()

    pool.close()
    pool.join()
    pbar.close()

    if unmapped_proteins:
        logging.info(f"Unable to process hits for protein(s) with following "
                     f"accession(s): \n {unmapped_proteins}\n"
                     f"Either the mapping to gene ID or the mapping of "
                     f"hit accession to respective taxon failed.")


def taxonomic_assignment(cfg, gff_df, assignments_df, tax_assignment_path, target_taxon, missing_taxids, TAX_DB):
    """
    Run functions to save taxonomic assignment to gene objects

    Args:
      tax_assignment_path:
      genes:
      prots:
      target_id:
      :param missing_taxids:
      :param TAX_DB:
    """

    process_hits(cfg, tax_assignment_path, assignments_df, target_taxon,
                 missing_taxids, TAX_DB)
    generalize_unclassified(assignments_df, missing_taxids, TAX_DB)


###############################################################################
########################### PROTEIN FASTA FILE ################################
###############################################################################


def subset_protein_fasta(proteins_path, prot_list, path_out, mode):
    """
    use the list of IDs to write to the tmp prot fasta only those
    proteins whose ID appears (not) in the list

    Args:
      proteins_path:
      prot_list:
      path_out:
      mode:

    Returns:

    """

    # mode: include or exclude list
    if mode == 'include':
        mode = True
    else:
        mode = False

    with open(proteins_path, 'r') as file_prot:
        with open(path_out, 'w') as file_out:
            write_bool = False
            for line in file_prot:
                if line.startswith('>'):
                    transcript = line.split()[0][1:]
                    # gene = line.split()[0].split('|')[1].lstrip(">")

                    if ((transcript in prot_list) and mode) or (
                            (transcript not in prot_list) and not mode):
                        write_bool = True
                        file_out.write(line)
                    else:
                        write_bool = False
                else:
                    if write_bool:
                        file_out.write(line)


def subset_prots_longest_cds(gff_df, proteins_path, path_out):
    """
    the gene-protein ID matching table, generated with gffread is used for this table
    (it contains the length of the CDS ) to infer the transcript with the longest cds for each gene
    this transcript is then written to a tmp fasta file

    Args:
      genes:
      proteins_path:
      path_out:

    Returns:

    """

    # when matching prot ID to gene ID it is already checked for the one with the longest CDS
    longest_transcripts = set(gff_df.loc[gff_df['type'] == 'gene', 'diamond_header'])

    logging.info(
        f"{len(longest_transcripts)} proteins written to fasta file for "
        f"taxonomic assignment (subsetting for longest CDS)")
    subset_protein_fasta(proteins_path, longest_transcripts, path_out,
                         "include")


def detect_appropriate_assignment(gene, quick_mode_match_id, missing_taxids, TAX_DB):


    tax_lineage = init_taxopyon(gene.taxon_assignmentID, missing_taxids,
                                TAX_DB).taxid_lineage
    # if quick_mode_match_id not in tax_lineage or quick_mode_match_id == tax_lineage[0]: # <- option 2: exclusive
    if quick_mode_match_id not in tax_lineage:  # <- option 1: inclusive
        # inclusive: if match id equals class, assignment must be at least same class as query
        # exclusive: if match id equals class, assignment must be at least one rank closer to query
        return gene.transcript_id
    else:
        return None




def filter_prots_quick_hits(assignments_df, quick_mode_match_id, prot_path_in,
                            prot_path_out, missing_taxids, TAX_DB):
    """
    write proteins to tmp protein FASTA which were not assigned in first
    quick search with DIAMOND

    Args:
      genes:
      quick_mode_match_id:
      prot_path_in:
      prot_path_out:

    Returns:

    """

    no_match = set(assignments_df.apply(
        lambda row: detect_appropriate_assignment(row,
                                                  quick_mode_match_id, missing_taxids, TAX_DB),
        axis=1))

    logging.info(f"{len(no_match)} proteins written to fasta file for 2nd DIAMOND run")
    subset_protein_fasta(prot_path_in, no_match, prot_path_out, "include")


###############################################################################
################################### MISC ######################################
###############################################################################


def diamond_inclusion_by_exclusion(include_id, missing_taxids, TAX_DB):
    """
    diamond can not exclude and include taxa simultaneously
    thus, this script identifies all taxIDs to exclude to simulate inclusion

    Args:
      include_id:

    Returns:
    :param TAX_DB:

    """

    exclude_list = []

    include_taxon = init_taxopyon(include_id, missing_taxids, TAX_DB)
    include_lineage = include_taxon.taxid_lineage
    for parent in include_lineage[1:]:
        childs = get_children_ids(TAX_DB, parent)
        exclude_list += [str(child) for child in childs if
                                       child not in include_lineage]

    return exclude_list


def taxonomy_summary(cfg, assignments_df, target_taxon, query_label, missing_taxids, TAX_DB):
    """
    computing a score to assess the quality of the assembly

    Args:
      cfg:
      assignments_df:

    Returns:
    :param TAX_DB:
    :param missing_taxids:

    """

    query_lineage = target_taxon.taxid_lineage
    total_assignments = 0
    outside_query_lineage = 0
    tax_sum_dict = {'superkingdom': {}, 'kingdom': {}, 'phylum': {},
                    'class': {}}  # , 'order': {}, 'family': {}, 'genus': {}, 'species': {}}
    unique_assignments = assignments_df.value_counts('taxon_assignmentID').to_frame(
        'count').reset_index()
    total_assignments = unique_assignments.shape[0]

    for index, row in unique_assignments.iterrows():
        g_taxon = init_taxopyon(row['taxon_assignmentID'], missing_taxids,
                                TAX_DB)
        if not g_taxon.taxid in query_lineage:
            outside_query_lineage += 1
        for rank in ['superkingdom', 'kingdom', 'phylum',
                     'class']:  # , 'order', 'family', 'genus', 'species']:
            rank_value = g_taxon.rank_name_dictionary.get(rank)
            if rank_value in tax_sum_dict[rank].keys():
                tax_sum_dict[rank][rank_value] = tax_sum_dict[rank].get(
                    rank_value) + row['count']
            else:
                tax_sum_dict[rank][rank_value] = row['count']

    with open(cfg.output_path + 'taxonomic_assignment/summary.txt',
              'w') as summary_file:
        summary_file.write(
            f'Query taxon: {cfg.taxon_name} (NCBI ID: {cfg.taxon_id})\n')
        summary_file.write(
            f'Assigend label for query: {query_label}\n\n')

        summary_file.write(f'Total number of distinct taxonomic '
                           f'assignments: {total_assignments}\n')
        summary_file.write(
            f'Taxonomic assignments not in target lineage: '
            f'{outside_query_lineage} '
            f'({round((outside_query_lineage/total_assignments)*100, 2)}%)\n\n')
        for rank, taxa in tax_sum_dict.items():
            summary_file.write(f'Following taxa for rank "{rank}" '
                               f'are present in the assembly ("taxon"\t"gene count"):\n')
            for taxon, count in taxa.items():
                if taxon != 'None' and taxon:
                    summary_file.write(f'{taxon}:\t{count}\n')


def str2float_or_int(string, digits):
    try:
        integer_number = int(string)
    except:
        integer_number = None
    try:
        float_number = round(float(string), digits)
    except:
        float_number = None

    if integer_number == None and float_number == None:
        return string
    elif integer_number != None:
        if float_number == integer_number:
            return integer_number
        else:
            return float_number
    else:
        return float_number


def write_output(cfg, assignments_df):
    """
    Write information in gene objects to gene_table_taxon_assignment.csv.

    Args:
      cfg:
      genes:
      header:

    Returns:

    """

    out_path = cfg.output_path + 'taxonomic_assignment/taxonomic_assignments.csv'
    assignments_df.to_csv(out_path, index=True, index_label='gene_id')

def assignment_lineages(assignments_df, missing_taxids, TAX_DB):

    rank_subset = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                   'family', 'genus', 'species', 'strain']

    for assignment in assignments_df['taxon_assignmentID'].dropna().unique():
        assigned_taxon = init_taxopyon(assignment, missing_taxids, TAX_DB)
        ranks_of_interest = {key: assigned_taxon.rank_name_dictionary.get(key)
                             for key in rank_subset}


###############################################################################
###############################################################################


def run_assignment(cfg, gff_df, pca_coordinates, TAX_DB):
    """

    Args:
      cfg:

    Returns:

    """

    missing_taxids = set()


    if cfg.assignment_mode == "quick":
        tax_assignment_path_1, tax_assignment_path_2 = cfg.tax_assignment_path
    else:
        tax_assignment_path = cfg.tax_assignment_path


    tmp_prot_path = cfg.output_path + "tmp/tmp.subset.protein.fasta"


    diamond_cmd = [cfg.diamond, 'blastp', '-p', str(cfg.threads), '-f', '6',
                   'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                   'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'sscinames',
                  '-b', '2.0', '--tmpdir', '/dev/shm', '-c', '1',
                  '--top', '10', '-q', tmp_prot_path, '-d', cfg.database_path]
    if cfg.diamond_sensitivity != 'default':
        diamond_cmd.append(f'--{cfg.diamond_sensitivity}')

    target_taxon = init_taxopyon(cfg.taxon_id, missing_taxids, TAX_DB)

    excluded_ids = []
    if cfg.target_exclude:
        excluded_ids.append(str(get_id_for_rank_of_species(target_taxon.taxid, cfg.exclusion_rank, missing_taxids, TAX_DB)))
    if cfg.add_exclude:
        if type(cfg.add_exclude) == list:
            for id in cfg.add_exclude:
                excluded_ids.append(str(id))
        elif ',' in cfg.add_exclude:
            for id in cfg.add_exclude.split(','):
                excluded_ids.append(str(id))
        else:
            excluded_ids.append(str(cfg.add_exclude))
    if cfg.target_exclude or cfg.add_exclude:
        target_exclude = ['--taxon-exclude'] + excluded_ids
    else:
        target_exclude = []


    if type(cfg.num_groups_plot) == int:
        num_groups_plot = cfg.num_groups_plot
    elif cfg.num_groups_plot.isdigit():
        num_groups_plot = int(cfg.num_groups_plot)
    elif cfg.num_groups_plot == 'all':
        # no merging of taxonomic assignments desired
        num_groups_plot = False
    else:
        logging.warning('No valid option for "num_groups_plot"')
        sys.exit()

    perform_diamond = cfg.compute_tax_assignment
    # let cfg.force not affect this
    # cfg.compute_tax_assignment is true when either file does not exist or
    # user specified to run it

    # init empty dataframe for taxonomic assignments
    num_genes = gff_df.loc[gff_df['type'] == 'gene'].shape[0]
    assignments_dict = {
        'diamond_header': gff_df.loc[gff_df['type'] == 'gene', 'diamond_header'].to_list(),
        'fasta_header': gff_df.loc[gff_df['type'] == 'gene', 'fasta_header'].to_list(),
        'best_hit': [None for i in range(num_genes)],
        'best_hitID': [None for i in range(num_genes)],
        'bh_pident': [None for i in range(num_genes)],
        'bh_evalue': [None for i in range(num_genes)],
        'bh_bitscore': [None for i in range(num_genes)],
        'lca': [None for i in range(num_genes)],
        'lcaID': [None for i in range(num_genes)],
        'refined_lca': [None for i in range(num_genes)],
        'refined_lcaID': [None for i in range(num_genes)],
        'taxon_assignment': [None for i in range(num_genes)],
        'taxon_assignmentID': [None for i in range(num_genes)],
        'is_target': [0 for i in range(num_genes)]
    }
    assignments_df = pd.DataFrame(assignments_dict,
        index = gff_df.loc[gff_df['type'] == 'gene'].index)

    if perform_diamond:
        subset_prots_longest_cds(gff_df, cfg.proteins_path, tmp_prot_path)
    if cfg.assignment_mode == 'quick':
        perform_quick_search_1(cfg, perform_diamond, diamond_cmd,
                               cfg.quick_mode_search_rank, target_exclude,
                               cfg.taxon_id, missing_taxids, TAX_DB)
        taxonomic_assignment(cfg, gff_df, assignments_df, tax_assignment_path_1,
                             target_taxon, missing_taxids, TAX_DB)
        perform_quick_search_2(cfg, perform_diamond, diamond_cmd,
                               cfg.quick_mode_match_rank, assignments_df,
                               tmp_prot_path, target_exclude, target_taxon,
                               missing_taxids, TAX_DB)
        taxonomic_assignment(cfg, gff_df, assignments_df, tax_assignment_path_2,
                             target_taxon, missing_taxids, TAX_DB)
    elif cfg.assignment_mode == 'exhaustive':
        perform_exhaustive_search(cfg, perform_diamond, diamond_cmd, target_exclude)
        taxonomic_assignment(cfg, gff_df, assignments_df, tax_assignment_path,
                             target_taxon, missing_taxids, TAX_DB)
    else:
        logging.error('Assignment mode not one of quick or exhaustive')
        sys.exit()

    assignments_df['plot_label'] = None
    assignments_df['plot_labelID'] = None
    assignments_df = merge_labels(assignments_df, target_taxon, cfg.merging_labels,
                 missing_taxids, TAX_DB)
    if num_groups_plot:
        num_labels = len(assignments_df['plot_labelID'].unique())
        assignments_df = iterative_merging_top(assignments_df, num_groups_plot - num_labels,
                              missing_taxids, TAX_DB)
    assignments_df = set_unassigned_labels(assignments_df)

    assignments_df = assignments_df.astype({'plot_labelID': 'Int64'})

    query_label = ident_query_label(assignments_df, target_taxon, TAX_DB)

    if missing_taxids:
        logging.info(
            "The following Taxon ID(s) could not be found in the NCBI (skipped in taxonomic assignment):")
        logging.info(missing_taxids)

    pathlib.Path(cfg.output_path + 'taxonomic_assignment/').mkdir(parents=True,
                                                                  exist_ok=True)

    #assignment_lineages(assignments_df)
    taxonomy_summary(cfg, assignments_df, target_taxon, query_label,
                     missing_taxids, TAX_DB)
    return assignments_df, query_label


def main():
    """ """
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    run_assignment(gff_df, cfg)


if __name__ == '__main__':
    main()
