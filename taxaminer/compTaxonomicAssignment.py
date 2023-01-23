#!/usr/bin/env python

"""Perform taxonomic assignment

Expects processed config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

from . import checkData

import taxopy
import sys
import csv
import pathlib
import subprocess
import logging
import json


class Gene:
    """ """

    def __init__(self, line_array, header_index, cov):
        self.g_name = line_array[header_index.get('g_name')]
        self.c_name = line_array[header_index.get('c_name')]
        self.c_num_of_genes = line_array[header_index.get('c_num_of_genes')]
        self.c_len = line_array[header_index.get('c_len')]
        self.c_pct_assemby_len = line_array[
            header_index.get('c_pct_assemby_len')]
        self.c_genelenm = line_array[header_index.get('c_genelenm')]
        self.c_genelensd = line_array[header_index.get('c_genelensd')]

        if cov:
            self.c_cov = {item_index: line_array[index] for item_index, index in
                          header_index.get('c_cov').items()}
            self.c_covsd = {item_index: line_array[index] for item_index, index
                            in header_index.get('c_covsd').items()}
            self.c_covdev = {item_index: line_array[index] for item_index, index
                             in header_index.get('c_covdev').items()}
            self.c_genecovm = {item_index: line_array[index] for
                               item_index, index in
                               header_index.get('c_genecovm').items()}
            self.c_genecovsd = {item_index: line_array[index] for
                                item_index, index in
                                header_index.get('c_genecovsd').items()}

        self.c_pearson_r = line_array[header_index.get('c_pearson_r')]
        self.c_pearson_p = line_array[header_index.get('c_pearson_p')]
        self.c_gc_cont = line_array[header_index.get('c_gc_cont')]
        self.c_gcdev = line_array[header_index.get('c_gcdev')]

        self.g_len = line_array[header_index.get('g_len')]
        self.g_lendev_c = line_array[header_index.get('g_lendev_c')]
        self.g_lendev_o = line_array[header_index.get('g_lendev_o')]
        self.g_abspos = line_array[header_index.get('g_abspos')]
        self.g_terminal = line_array[header_index.get('g_terminal')]
        self.g_single = line_array[header_index.get('g_single')]

        if cov:
            self.g_cov = {item_index: line_array[index] for item_index, index in
                          header_index.get('g_cov').items()}
            self.g_covsd = {item_index: line_array[index] for item_index, index
                            in header_index.get('g_covsd').items()}
            self.g_covdev_c = {item_index: line_array[index] for
                               item_index, index in
                               header_index.get('g_covdev_c').items()}
            self.g_covdev_o = {item_index: line_array[index] for
                               item_index, index in
                               header_index.get('g_covdev_o').items()}

        self.g_pearson_r_o = line_array[header_index.get('g_pearson_r_o')]
        self.g_pearson_p_o = line_array[header_index.get('g_pearson_p_o')]
        self.g_pearson_r_c = line_array[header_index.get('g_pearson_r_c')]
        self.g_pearson_p_c = line_array[header_index.get('g_pearson_p_c')]
        self.g_gc_cont = line_array[header_index.get('g_gc_cont')]
        self.g_gcdev_c = line_array[header_index.get('g_gcdev_c')]
        self.g_gcdev_o = line_array[header_index.get('g_gcdev_o')]

        self.coords = None

        self.fasta_header = None  # ID of transcript with longest CDS
        self.lencds = 0  # length of self.fasta_header
        self.prot_seq = None

        self.tax_assignments = []  # all taxonomic hits for the gene
        self.best_hit = None
        self.best_hitID = None
        self.bh_evalue = None
        self.bh_pident = None
        self.lca = None
        self.lcaID = None
        self.corrected_lca = None
        self.corrected_lcaID = None
        self.taxon_assignment = None
        self.taxon_assignmentID = None

        self.plot_label = None
        self.plot_labelID = None


###############################################################################
############################ HELPER FUNCTIONS #################################
###############################################################################

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
    for prefix in ['gene:', 'gene-', 'transcript:', 'transcript-',
                   'rna:', 'rna-', 'cds:', 'cds-']:
        id = remove_prefix(id, prefix)
    return id


def get_id_for_rank_of_species(queryID, rank):
    """
    Returns ID of taxon at given rank for given taxon ID.

    Args:
      queryID(int): NCBI ID of taxon for wich to find ID at given rank
      rank(str): rank of which to get the NCBI ID of

    Returns:
       (int) NCBI ID of taxon at the given rank for given taxon ID
     """
    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_ranks = query_taxon.rank_name_dictionary
    rank_of_interest = query_ranks.get(rank)
    # get ID for the specific class (function taxopy.taxid_from_name() might get multiple IDs because of homonyms)
    if rank_of_interest:
        rank_index = query_taxon.name_lineage.index(rank_of_interest)
        rank_id = query_taxon.taxid_lineage[rank_index]
        return rank_id
    else:
        return None


def get_children_ids(parent):
    """
    Returns list of taxon IDs which are children of given taxon

    Args:
      parent (int): NCBI taxon ID to get IDs of children for

    Returns:
        (list) list of taxon IDs which are children of given taxon
    """
    # taxid2parent: dict with key = child, value = parent
    return [key for key in TAX_DB.taxid2parent.keys() if
            TAX_DB.taxid2parent[key] == parent]


def get_gff_attribute(attr_list, attr):
    """
    Return value for GFF attribute of given list of attributes

     Args:
      attr_list:
      attr:

    Returns:

    """
    tmp = attr_list.split(attr + '=')[1]
    if len(tmp.split(';')) > 1:
        value = strip_ID(tmp.split(';')[0])
    else:
        value = strip_ID(tmp)
    return value


###############################################################################
######################### PLOT LABEL COMPUTATION ##############################
###############################################################################


def init_label_taxa_dict(genes):
    """
    initializes dictionary which stores which taxonomic assignments are summarized in each label
    every taxonomic assignment is assigned to 'root' in the beginning
    # also set initial plot label and ID (to root) for each gene

    Args:
      genes:

    Returns:

    """

    labelID_taxa = {1: set()}  # taxa that have this label; label : set(taxa)
    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            gene.plot_label = 'Unassigned'
            gene.plot_labelID = 'NA'
            continue
        labelID_taxa[1].add(gene.taxon_assignmentID)
    return labelID_taxa


def get_children_for_label(taxa, labelID):
    """
    gets all descendants of label which are in lineage of one of the taxa the label represents
    also counts how many of these descendants are already in list of current labels

    Args:
      taxa:
      labelID:

    Returns:

    """
    children = {}
    for taxID in taxa:
        taxon = taxopy.Taxon(taxID, TAX_DB)
        child = taxon.taxid_lineage[max(taxon.taxid_lineage.index(labelID) - 1,
                                        0)]  # max 0 if taxon is end of lineage
        if child in children:
            children[child].add(taxID)
        else:
            children[child] = set([taxID])

    return children


def merge_assignments_at_rank(genes, rank):
    """
    Merge taxonomic assignments at given rank.
    Args:
      genes(dict): dictionary of genes {g_name: gene object}
      rank(str): rank at which to merge all taxonomic assignments at
    """
    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            gene.plot_label = 'Unassigned'
            gene.plot_labelID = 'NA'
            continue
        taxon = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB)
        ranks = taxon.rank_name_dictionary
        merger_rank = ranks.get(rank)
        if merger_rank:
            gene.plot_label = merger_rank
            gene.plot_labelID = get_id_for_rank_of_species(
                gene.taxon_assignmentID, rank)
        else:
            gene.plot_label = gene.taxon_assignment
            gene.plot_labelID = gene.taxon_assignmentID


def merge_assignments_at_id(genes, ids):
    """
    Merge taxonomic assignments at given NCBI taxon IDs.

    Args:
      genes(dict): dictionary of genes {g_name: gene object}
      ids(int,list): list or int of taxon ID(s) at which to merge
    """

    if type(ids) != list:
        ids = [ids]

    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            gene.plot_label = 'Unassigned'
            gene.plot_labelID = 'NA'
            continue

        taxon = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB)
        lineage = taxon.taxid_lineage
        for id in ids:
            if id in lineage:
                id_merger_taxon = taxopy.Taxon(id, TAX_DB)
                gene.plot_label = id_merger_taxon.name
                gene.plot_labelID = id_merger_taxon.taxid
                lineage = id_merger_taxon.taxid_lineage


def iterative_merging_top(genes, num_groups_plot):
    """
    identify at which rank in query lineage to merge best

    Args:
      genes:
      num_groups_plot:

    Returns:

    """

    labelID_taxa = init_label_taxa_dict(
        genes)  # taxa that have this label; label : set(taxa)

    update_dict = {'nonempty-dummy'}
    while update_dict:
        update_dict = {}
        del_dict_keys = set()
        remaining_labels_num = num_groups_plot - len(labelID_taxa)
        for labelID, taxa in labelID_taxa.items():
            if len(taxa) >= 2:
                lca_candidates = []
                for taxon in taxa:
                    lca_candidates.append(taxopy.Taxon(taxon, TAX_DB))
                lca = taxopy.find_lca(lca_candidates, TAX_DB)
                if lca.taxid != labelID:  # if it is possible to find a more specific label than the current one with the lca
                    del_dict_keys.add(labelID)
                    update_dict[lca.taxid] = taxa
                else:
                    # get all children for label that are in list taxa
                    children = get_children_for_label(taxa, labelID)
                    existing_labels = len(children.keys()) - len(
                        list(set(set(children.keys()) - labelID_taxa.keys())))

                    # check if number of new labels (= num of children - labels that are already in lablels)
                    # fits in remainung number of possible new labels
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
            taxid2label[min(taxa)] = min(taxa)
        else:
            for taxon in taxa:
                taxid2label[taxon] = labelID

    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            continue
        else:
            gene.plot_label = TAX_DB.taxid2name[
                taxid2label.get(gene.taxon_assignmentID)]
            gene.plot_labelID = taxid2label.get(gene.taxon_assignmentID)


def filter_labeled_genes(genes):
    """Filter genes which already have a label and return number of labels.

    Args:
      genes(dict): dictionary of genes {g_name: gene object}

    Returns:
      dict: dictionary of genes without label {g_name: gene object}
      int: number of distinct labels in the set of labeled genes
    """
    unlabeled_genes = {}
    labels = set()

    for g_name, gene in genes.items():
        if not gene.plot_label:
            unlabeled_genes[g_name] = gene
        else:
            labels.add(gene.plot_labelID)

    num_labels = len(labels)
    return unlabeled_genes, num_labels


def merge_labels(genes, queryID, merging_labels):
    """
    Merge plot labels based on user input.

    Args:
      genes(dict): dictionary of genes {g_name: gene object}
      queryID(int): NCBI taxon ID for query species
      merging_labels(list,int): user input for merging labels
        1. taxonomic rank with suffix '-all'
            -> merge all taxonomic assignment at that rank
        2. taxonomic rank
            -> merge query lineage at that rank
        3. list of NCBI taxon IDs
            -> merge at each of the taxa provided
    """

    merger_ids = []
    all_merger_rank = None
    query_merger_id = None

    if not merging_labels:
        merging_labels = []
    elif type(merging_labels) != list:
        merging_labels = [merging_labels]

    for merger_label in merging_labels:

        if not merger_label:
            continue

        if type(merger_label) == str:
            # option 1: string end with '-all' -> every taxonomic assignment is merged at that rank
            # option 2: string does not end with '-all' -> only taxonomic assignments with same rank <merger_label> as the query species are merged
            if merger_label.endswith('-all'):
                all_merger_rank = merger_label[:-len('-all')]
                break
            else:
                # get ncbiID for rank <merger_label> of query
                query_merger_id = get_id_for_rank_of_species(queryID,
                                                             merger_label)
                break
        else:
            if type(merger_label) == int or merger_label.isdigit():
                merger_ids.append(
                    int(merger_label))  # ncbiID of <merger_label> to merge at

    # if there is an -all merger-label; this overwrites everything
    # else; merge at query_merger_id and merger_ids
    if all_merger_rank:
        merge_assignments_at_rank(genes, all_merger_rank)
    elif query_merger_id:
        merge_assignments_at_id(genes, query_merger_id)
    elif merger_ids:
        merge_assignments_at_id(genes, merger_ids)


def set_unassigned_labels(genes):
    """
    Set missing plot labels to taxonomic assignment or 'Unassigned'.

    Args:
      genes:
    """
    for g_name, gene in genes.items():
        if not gene.plot_label:
            if gene.taxon_assignmentID == 'NA':
                gene.plot_label = 'Unassigned'
                gene.plot_labelID = 'NA'
            else:
                gene.plot_label = gene.taxon_assignment
                gene.plot_labelID = gene.taxon_assignmentID


def ident_query_label(genes, queryID):
    """
    Identify which label represents the query species.

    Args:
      genes:
      queryID:

    Returns:

    """

    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_lineage = query_taxon.taxid_lineage

    label_candidate = (None, len(query_lineage))

    for g_name, gene in genes.items():
        if gene.plot_labelID in query_lineage:
            if query_lineage.index(gene.plot_labelID) < label_candidate[1]:
                label_candidate = (
                gene.plot_label, query_lineage.index(gene.plot_labelID))

    return label_candidate[0]


###############################################################################
########################## TAXONOMIC ASSIGNMENT ###############################
###############################################################################

def get_subset_id(queryID, subset_marker):
    """
    Return input if input.isdigit() or return ID for rank of query.

    Args:
      queryID:
      subset_marker:

    Returns:

    """
    if type(subset_marker) == int:
        return subset_marker
    elif subset_marker.isdigit():
        return int(subset_marker)
    else:
        return get_id_for_rank_of_species(queryID, subset_marker)


def run_diamond(diamond_cmd):
    """
    perform DIAMOND run and catch errors

    Args:
      diamond_cmd:

    Returns:

    """
    dmdnOut = subprocess.run([diamond_cmd], shell=True, capture_output=True)
    if dmdnOut.returncode != 0:
        logging.error(f'Error running DIAMOND '
                      f'with following command:\n{diamond_cmd}')
        if '\nOpening the database... No such file or directory\n' in dmdnOut.stderr.decode():
            logging.error('Database does not exists. Please check your input.')
        # TODO: check for output when putting empty protein file
        elif 'chuighdseck' in dmdnOut.stderr.decode():
            logging.error('Empty protein FASTA file. Please check your input.')
        else:
            logging.error('DIAMOND Error message:\n' + dmdnOut.stderr.decode())
        sys.exit()
    else:
        logging.info('DIAMOND ouput:\n')
        for line in dmdnOut.stderr.decode().split('\n'):
            if line.startswith('Scoring parameters:'):
                logging.debug(line)
            elif line.startswith('Temporary directory:'):
                logging.debug(line)
            elif line.startswith(
                    'Percentage range of top alignment score to report hits:'):
                logging.debug(line)
            elif line.startswith('Reference ='):
                logging.debug(line)
            elif line.startswith('Sequences ='):
                logging.debug(line)
            elif line.startswith('Letters ='):
                logging.debug(line)
            elif line.startswith('Block size ='):
                logging.debug(line)
            elif line.startswith('Total time ='):
                logging.info(f'DIAMOND runtime: {line.split(" = ")[1]}')
            elif ' pairwise alignments, ' in line:
                logging.info(line)
            elif ' queries aligned.' in line:
                logging.info(line)


def perform_quick_search_1(perform_diamond, diamond_cmd, tax_assignment_path_1,
                           quick_mode_search_rank, taxon_exclude, queryID):
    """
    Perform first DIAMOND run for quick assignment mode.

    Args:
      perform_diamond:
      diamond_cmd:
      tax_assignment_path_1:
      quick_mode_search_rank:
      taxon_exclude:
      queryID:
    """
    logging.info("Quick mode for taxonomic assignment selected")

    if perform_diamond:
        # quick_mode_search_rank can be either taxon ID or rank
        quick_mode_search_id = get_subset_id(queryID, quick_mode_search_rank)

        if taxon_exclude:
            q1_exclude = taxon_exclude.rstrip('"') + ',' + ','.join(
                diamond_inclusion_by_exclusion(quick_mode_search_id))
        else:
            q1_exclude = f' --taxon-exclude ' \
                         f'"{",".join(diamond_inclusion_by_exclusion(quick_mode_search_id))}"'
        diamond_cmd_1 = f'{diamond_cmd} -o "{tax_assignment_path_1}" {q1_exclude}'
        ###taxonlist###  diamond_cmd_1 = diamond_cmd + diamond_o1 + ' --taxonlist "' + str(quick_mode_search_id) + '" --taxon-exclude "' + str(queryID) + '"'

        logging.debug(f"Diamond command for "
                      f"quick search round 1:\n{diamond_cmd_1}")
        run_diamond(diamond_cmd_1)


def perform_quick_search_2(perform_diamond, diamond_cmd, tax_assignment_path_2,
                           quick_mode_match_rank, genes, tmp_prot_path,
                           taxon_exclude, queryID):
    """
    Perform second DIAMOND run for quick assignment mode.

    Args:
      perform_diamond:
      diamond_cmd:
      tax_assignment_path_2:
      quick_mode_match_rank:
      genes:
      tmp_prot_path:
      taxon_exclude:
      queryID:
    """

    if perform_diamond:
        quick_mode_match_id = get_subset_id(queryID, quick_mode_match_rank)

        filter_prots_quick_hits(genes, quick_mode_match_id, tmp_prot_path,
                                tmp_prot_path)

        diamond_cmd_2 = f'{diamond_cmd} -o "{tax_assignment_path_2}" {taxon_exclude}'
        logging.debug(f"Diamond command for quick search round 2:\n{diamond_cmd_2}")
        run_diamond(diamond_cmd_2)


def perform_exhaustive_search(perform_diamond, diamond_cmd,
                              tax_assignment_path, taxon_exclude):
    """
    Perform  DIAMOND run for exhaustive assignment mode.

    Args:
      perform_diamond:
      diamond_cmd:
      tax_assignment_path:
      taxon_exclude:
    """
    logging.info("Exhaustive mode for taxonomic assignment selected")
    if perform_diamond:
        diamond_cmd = f'{diamond_cmd} -o "{tax_assignment_path}" {taxon_exclude}'
        logging.debug(
            f"Diamond command for exhaustive search:\n{diamond_cmd}")
        run_diamond(diamond_cmd)


def assess_best_of_multi_hits(taxIDs, queryID):
    """
    for a list of taxa IDs, get the one where LCA with queryID is closest to query

    Args:
      taxIDs:
      queryID:

    Returns:

    """

    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_lineage = query_taxon.taxid_lineage

    bestID = (None, len(query_lineage) + 1)  # (id, index in lineage)
    for taxID in taxIDs:
        if taxID in TAX_DB.taxid2name.keys():  # avoid crashes because of ncbiIDs not being represented in TAX_DB
            taxon = taxopy.Taxon(taxID, TAX_DB)
            lca = taxopy.find_lca([query_taxon, taxon], TAX_DB)
            lca_index = query_lineage.index(lca.taxid)
            if bestID[
                1] > lca_index:  # check if new best hit is closer to query in lineage as current best hit
                bestID = (taxon, lca_index)
        else:
            missing_taxids.add(taxID)

    return bestID[0]


def set_taxon_assignment(genes):
    """
    Assess taxonomic assignment for each gene.

    If there is a corrected_lca; taxonomic assignment = corrected_lca
    else; taxonomic assignment = lca

    Args:
      genes:
    """
    for g_name, gene in genes.items():
        if gene.corrected_lca != 'NA':
            gene.taxon_assignment = gene.corrected_lca
            gene.taxon_assignmentID = gene.corrected_lcaID
        else:
            gene.taxon_assignment = gene.lca
            gene.taxon_assignmentID = gene.lcaID


def calc_corrected_lca(genes, queryID):
    """
    Calculate the corrected LCA for each gene based on best hit and query.

    Args:
      genes:
      queryID:
    """
    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_lineage = query_taxon.taxid_lineage

    for g_name, gene in genes.items():
        if gene.lcaID in query_lineage and gene.best_hitID != 'NA':
            if gene.best_hitID in TAX_DB.taxid2name.keys():  # avoid crashes because of ncbiIDs not being represented in TAX_DB
                bh_taxon = taxopy.Taxon(gene.best_hitID, TAX_DB)
                corrected_lca = taxopy.find_lca([query_taxon, bh_taxon], TAX_DB)
                gene.corrected_lcaID = corrected_lca.taxid
                gene.corrected_lca = corrected_lca.name
            else:
                missing_taxids.add(gene.best_hitID)
                gene.corrected_lcaID = 'NA'
                gene.corrected_lca = 'NA'
        else:
            gene.corrected_lcaID = 'NA'
            gene.corrected_lca = 'NA'


def calc_lca(genes):
    """
    Calculate the LCA for each gene based on taxonomic hits.

    Args:
      genes:

    Returns:

    """

    for g_name, gene in genes.items():
        lca_query = []
        for tax_assignment in gene.tax_assignments:
            if tax_assignment[-2]:
                tax_id = int(tax_assignment[-2].split(';')[0])
                if tax_id in TAX_DB.taxid2name.keys():  # avoid crashes because of ncbiIDs not being represented in TAX_DB
                    taxon = taxopy.Taxon(tax_id, TAX_DB)
                    lca_query.append(taxon)
                else:
                    missing_taxids.add(tax_id)
        if len(lca_query) > 1:
            lca = taxopy.find_lca(lca_query, TAX_DB)
            gene.lca = lca.name
            gene.lcaID = int(lca.taxid)
        elif not lca_query:
            gene.lca = 'NA'
            gene.lcaID = 'NA'
            pass
        else:  # if only one hit for taxon
            gene.lca = lca_query[0].name
            gene.lcaID = int(lca_query[0].taxid)


def read_tax_assignments(tax_assignment_path, prots, queryID):
    """
    Read DIAMOND output file and assign hits to genes.

    Args:
      tax_assignment_path:
      prots:
      queryID:

    Returns:

    """

    with open(tax_assignment_path, 'r') as tax_assignment:
        # cols: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore taxid taxname
        for line in tax_assignment:
            # first is "best hit" in terms of alignment score
            spline = line.strip().split('\t')
            gene = prots.get(spline[0])
            closest_hit = None
            if gene:
                if not gene.best_hit or gene.best_hit == 'NA':  # if first best hit had no matching to NCBI ID
                    if spline[-2]:
                        if ';' in spline[-2]:
                            hit_ids = [int(id) for id in spline[-2].split(';')]
                            closest_hit = assess_best_of_multi_hits(hit_ids,
                                                                    queryID)
                            if closest_hit:
                                gene.best_hit = closest_hit.name
                                gene.best_hitID = closest_hit.taxid
                            else:
                                gene.best_hit = 'NA'
                                gene.best_hitID = 'NA'
                        else:
                            gene.best_hit = spline[-1]
                            gene.best_hitID = int(spline[-2])
                        gene.bh_evalue = spline[10]
                        gene.bh_pident = spline[2]
                    else:  # ncbiID for match could not be found
                        gene.best_hit = 'NA'
                        gene.best_hitID = 'NA'
                        gene.bh_evalue = 'NA'
                        gene.bh_pident = 'NA'

                if closest_hit:
                    gene.tax_assignments.append(
                        spline[1:-2] + [str(closest_hit.taxid),
                                        closest_hit.name])
                elif ';' in spline[-2]:
                    hit_ids = [int(id) for id in spline[-2].split(';')]
                    closest_hit = assess_best_of_multi_hits(hit_ids, queryID)
                    if closest_hit:
                        gene.tax_assignments.append(
                            spline[1:-2] + [str(closest_hit.taxid),
                                            closest_hit.name])
                    else:
                        gene.tax_assignments.append(spline[1:])
                else:
                    gene.tax_assignments.append(spline[1:])


def taxonomic_assignment(tax_assignment_path, genes, prots, queryID):
    """
    Run functions to save taxonomic assignment to gene objects

    Args:
      tax_assignment_path:
      genes:
      prots:
      queryID:
    """
    read_tax_assignments(tax_assignment_path, prots, queryID)
    calc_lca(genes)
    calc_corrected_lca(genes, queryID)
    set_taxon_assignment(genes)


def reset_tax_assignment(gene):
    """
    Delete all taxonomic assignment info for given gene

    Args:
      gene:
    """
    gene.tax_assignments = []
    gene.best_hit = None
    gene.best_hitID = None
    gene.lca = None
    gene.lcaID = None
    gene.corrected_lca = None
    gene.corrected_lcaID = None
    gene.taxon_assignment = None
    gene.taxon_assignmentID = None


###############################################################################
############################# DATA PREPARATION ################################
###############################################################################


def get_child_parent(attrs, child_attr, parent_attr):
    """
    Use GFF attributes and retrieve child-parent pairs, given certain attributes.

    Args:
      attrs:
      child_parent_dict:

    Returns:

    """
    if child_attr + "=" in attrs and parent_attr + "=" in attrs:
        childID = get_gff_attribute(attrs, child_attr)
        parentID = get_gff_attribute(attrs, parent_attr)
    else:
        return None, None
    return childID, parentID


def add_to_dict(dict, key, value):
    """
    add key value pair to dict
    and replace entry of key if value for key equals key

    Args:
      dict:
      key:
      value:
    """
    if not dict.get(key):
        dict[key] = value
    elif dict.get(key) == key:
        dict[key] = value


def traverse_dict(parent, child_parent_dict):
    """
    traverse trough child_parent_dict to find final parent

    Args:
      parent:
      child_parent_dict:

    Returns:

    """
    while parent in child_parent_dict.keys():
        if parent == child_parent_dict.get(parent):
            break
        else:
            parent = child_parent_dict.get(parent)
    return parent


def assign_prot2gene(prots, gene, id, length):
    """
    assign protein ID to gene ID whilst checking for CDS length

    Args:
      prots:
      gene:
      id:
      length:
    """
    if gene.fasta_header:
        if length > gene.lencds:
            gene.fasta_header = id
            gene.lencds = length
            prots[gene.fasta_header] = gene
    else:
        gene.fasta_header = id
        gene.lencds = length
        prots[gene.fasta_header] = gene


def prot_gene_matching(output_path, gff_path, genes, cfg, omitted_genes):
    """
    Match headers in protein FASTA to gene ID.
    Parses the GFF using the given gff parsing rule

    Args:
      output_path:
      gff_path:
      genes:
      cfg:

    Returns:

    """
    prots = {}

    child_parent_dict = {}  # key = child, value = parent
    with open(gff_path, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('#'):
                spline = line.strip().split('\t')

                # read lines where fasta headers are contained and/or child-parent-relations are to be found
                if cfg.gff_source == "default" or spline[1] == cfg.gff_source:
                    if spline[2] in cfg.gff_fasta_header_type:
                        if cfg.gff_connection == "parent/child":
                            # when attribute for child-parent-matching unequals attribute where header can be found,
                            # match attribute of fasta header and child of parent/child matching together
                            if cfg.gff_child_attr not in cfg.gff_fasta_header_attr:
                                childID, parentID = get_child_parent(spline[8],
                                                                     cfg.gff_fasta_header_attr,
                                                                     cfg.gff_child_attr)
                                if childID:
                                    add_to_dict(child_parent_dict,
                                                strip_ID(childID),
                                                strip_ID(parentID))
                            # add the parent child matching if fasta header line (also) contains p/c matching info
                            if spline[2] in cfg.gff_parent_child_types:
                                childID, parentID = get_child_parent(spline[8],
                                                                     cfg.gff_child_attr,
                                                                     cfg.gff_parent_attr)
                                if childID:
                                    add_to_dict(child_parent_dict,
                                                strip_ID(childID),
                                                strip_ID(parentID))
                        elif cfg.gff_connection == "inline":
                            headerID = get_gff_attribute(spline[8],
                                                         cfg.gff_fasta_header_attr)
                            geneID = get_gff_attribute(spline[8],
                                                       cfg.gff_gene_attr)
                            add_to_dict(child_parent_dict, strip_ID(headerID),
                                        strip_ID(geneID))
                    elif cfg.gff_connection == "parent/child" and spline[
                        2] in cfg.gff_parent_child_types:
                        childID, parentID = get_child_parent(spline[8],
                                                             cfg.gff_child_attr,
                                                             cfg.gff_parent_attr)
                        if childID:
                            add_to_dict(child_parent_dict, strip_ID(childID),
                                        strip_ID(parentID))

            elif "#FASTA" in line:  # if FASTA block has been reached
                break

    prot_index_path = output_path + 'tmp/tmp.proteins.fa.fai'

    # print(child_parent_dict)

    patternmatch = []
    unmatched = []  # list of proteins where gene could not be matched

    with open(prot_index_path, 'r') as prot_index:
        for line in prot_index:
            id, length = line.split()[0], ((int(line.split()[1]) * 3) + 3)
            parent = child_parent_dict.get(strip_ID(id))

            if parent:  # id in proteins fasta index matches the IDs in the GFF
                parent = traverse_dict(parent, child_parent_dict)
            elif child_parent_dict.get(
                    strip_ID(id.split('|')[0])):  # header might contain pipes
                parent = child_parent_dict.get(id.split('|')[0])
                parent = traverse_dict(parent, child_parent_dict)
            elif strip_ID(id) in genes.keys():
                # fasta header is gene ID
                parent = strip_ID(id)
            elif strip_ID(id).split('|')[0] in genes.keys():
                parent = strip_ID(id).split('|')[0]

            if not parent:
                patternmatch.append((id, length))
            else:
                gene = genes.get(parent)
                # put the ID as the fasta_header, as this will be what is used in DIAMOND results
                if gene:
                    assign_prot2gene(prots, gene, id, length)

    for id, length in patternmatch:
        # look if any child ID is contained in the id
        parent = None
        for child in child_parent_dict.keys():
            if child in id.split('|'):
                parent = child_parent_dict.get(child)
                parent = traverse_dict(parent, child_parent_dict)
        if not parent:  # still no success try to see if you can find gene ID in id
            for g_name in genes.keys():
                if g_name in id.split('|'):
                    parent = g_name

        if not parent and (not id in omitted_genes):
            unmatched.append(id)
        else:
            gene = genes.get(parent)
            # put the ID as the fasta_header, as this will be what is used in DIAMOND results
            # gene has to exist and protein ID is either not assigned yet or new ID is same as old one
            # (prevent overwriting of errorneous IDs)
            if gene and (not gene.fasta_header or gene.fasta_header == id):
                assign_prot2gene(prots, gene, id, length)

    if len(prots) == 0:
        logging.error(
            "Proteins were unable to be matched with gene IDs via GFF. Please "
            "check that FASTA headers fit with IDs in GFF (string before first "
            "whitespace or pipe must match ID of mRNA or CDS of corresponding gene)")
        sys.exit()
    else:
        if unmatched:
            logging.warning(f"Protein(s) with following FASTA header could not "
                            f"be matched to a gene ID:\n{unmatched}")
    return prots


def parse_header(header):
    """
    Returns dictionary mapping column names to index in the header of gene_table_coords.

    Args:
      header:

    Returns:

    """

    indexing = {}
    for index, item in enumerate(header.split(',')):
        raw_item = item.strip().rstrip('_0123456789')
        if raw_item in indexing.keys():  # for coverage variables
            item_index = item.strip().split('_')[-1].split('.')[-1]
            indexing[raw_item][item_index] = index
        elif raw_item == 'Dim.' or 'cov' in raw_item:
            item_index = item.strip().split('_')[-1].split('.')[-1]
            indexing[raw_item] = {item_index: index}
        else:
            indexing[raw_item] = index
    return indexing


def read_raw_genes(cfg):
    """
    Read gene_table_coords.csv and save information to Gene objects.

    Args:
      output_path:

    Returns:

    """

    input_path = cfg.output_path + 'gene_info/raw_gene_table.csv'

    with open(input_path, 'r') as input_file:
        header = next(input_file).strip()
        header_index = parse_header(header)

        genes = {}
        for line in input_file:
            gene = Gene(line.strip().split(','), header_index,
                        cfg.include_coverage)
            genes[gene.g_name] = gene

    return genes, header


def assign_pca_coords(genes, pca_coordinates, header):
    omitted_genes = []
    for g_name, gene in genes.items():
        gene.coords = {i: pca_coordinates.loc[g_name, f'PC {i}'] for i in
                       range(1, len(pca_coordinates.columns) + 1)}

        # filter genes without PCA coordinates
        if list(gene.coords.values())[0] == "NA":
            omitted_genes.append(gene.g_name)
            del genes[g_name]

    header = header + ',' + ','.join(
        f'PC_{i}' for i in range(1, len(pca_coordinates.columns) + 1))

    return header, omitted_genes


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
        lines = file_prot.readlines()
    with open(path_out, 'w') as file_out:
        write_bool = False
        for line in lines:
            if line.startswith('>'):
                transcript = line.split()[0].lstrip(">")
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


def subset_prots_longest_cds(genes, proteins_path, path_out):
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
    longest_transcripts = [gene.fasta_header for gene in genes.values() if
                           gene.fasta_header]

    logging.info(
        f"{len(longest_transcripts)} proteins written to fasta file for "
        f"taxonomic assignment (subsetting for longest CDS)")
    subset_protein_fasta(proteins_path, longest_transcripts, path_out,
                         "include")


def filter_prots_quick_hits(genes, quick_mode_match_id, prot_path_in,
                            prot_path_out):
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

    no_match = []

    for g_name, gene in genes.items():
        if gene.fasta_header:
            if gene.taxon_assignmentID != 'NA':

                tax_lineage = taxopy.Taxon(gene.taxon_assignmentID,
                                           TAX_DB).taxid_lineage
                # if quick_mode_match_id not in tax_lineage or quick_mode_match_id == tax_lineage[0]: # <- option 2: exclusive
                if quick_mode_match_id not in tax_lineage:  # <- option 1: inclusive
                    # inclusive: if match id equals class, assignment must be at least same class as query
                    # exclusive: if match id equals class, assignment must be at least one rank closer to query
                    no_match.append(gene.fasta_header)
                    reset_tax_assignment(gene)
            else:
                no_match.append(gene.fasta_header)
                reset_tax_assignment(gene)

    logging.info(f"{len(no_match)} proteins written to fasta file for 2nd DIAMOND run")
    subset_protein_fasta(prot_path_in, no_match, prot_path_out, "include")


###############################################################################
################################### MISC ######################################
###############################################################################


def diamond_inclusion_by_exclusion(include_id):
    """
    diamond can not exclude and include taxa simultaneously
    thus, this script identifies all taxIDs to exclude to simulate inclusion

    Args:
      include_id:

    Returns:

    """

    exclude_list = []

    include_taxon = taxopy.Taxon(include_id, TAX_DB)
    include_lineage = include_taxon.taxid_lineage
    for parent in include_lineage[1:]:
        childs = get_children_ids(parent)
        exclude_list = exclude_list + [str(child) for child in childs if
                                       child not in include_lineage]

    return exclude_list


def taxonomy_summary(cfg, genes):
    """
    computing a score to assess the quality of the assembly

    Args:
      cfg:
      genes:

    Returns:

    """

    query_taxon = taxopy.Taxon(cfg.taxon_id, TAX_DB)
    query_label = ident_query_label(genes, cfg.taxon_id)
    query_name = TAX_DB.taxid2name[cfg.taxon_id]
    query_lineage = query_taxon.taxid_lineage
    total_assignments = 0
    outside_query_lineage = 0
    tax_sum_dict = {'superkingdom': {}, 'kingdom': {}, 'phylum': {},
                    'class': {}}  # , 'order': {}, 'family': {}, 'genus': {}, 'species': {}}
    for gene in genes.values():
        if gene.taxon_assignmentID and gene.taxon_assignmentID != 'NA':
            total_assignments += 1
            g_taxon = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB)
            if not g_taxon.taxid in query_lineage:
                outside_query_lineage += 1
            for rank in ['superkingdom', 'kingdom', 'phylum',
                         'class']:  # , 'order', 'family', 'genus', 'species']:
                rank_value = g_taxon.rank_name_dictionary.get(rank)
                if rank_value in tax_sum_dict[rank].keys():
                    tax_sum_dict[rank][rank_value] = tax_sum_dict[rank].get(
                        rank_value) + 1
                else:
                    tax_sum_dict[rank][rank_value] = 1

    with open(cfg.output_path + 'taxonomic_assignment/summary.txt',
              'w') as summary_file:
        summary_file.write(
            f'Query taxon: {query_name} (NCBI ID: {cfg.taxon_id})\n')
        summary_file.write(
            f'Assigend label for query: {query_label}\n\n')

        summary_file.write(f'Total number of genes with taxonomic '
                           f'assignments: {total_assignments}\n')
        summary_file.write(
            f'Taxonomic assignments outside of query lineage: '
            f'{outside_query_lineage} '
            f'({round((outside_query_lineage/total_assignments)*100, 2)}%)\n\n')
        for rank, taxa in tax_sum_dict.items():
            summary_file.write(f'Following taxa for rank "{rank}" '
                               f'are present in the assembly:\n')
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


def write_output(cfg, genes, header):
    """
    Write information in gene objects to gene_table_taxon_assignment.csv.

    Args:
      cfg:
      genes:
      header:

    Returns:

    """

    out_path = cfg.output_path + 'taxonomic_assignment/gene_table_taxon_assignment.csv'

    cov_columns = ['c_cov_', 'c_covsd_', 'c_covdev_', 'c_genecovm_',
                   'c_genecovsd_', 'g_cov_', 'g_covsd_', 'g_covdev_c_',
                   'g_covdev_o_']

    with open(out_path, 'w') as out_file:

        csv_columns = header.strip().split(',') + ['fasta_header', 'lcaID',
                                                   'lca', \
                                                   'best_hitID', 'best_hit',
                                                   'bh_evalue', 'bh_pident', \
                                                   'corrected_lca',
                                                   'taxon_assignment',
                                                   'taxon_assignmentID',
                                                   'plot_label']
        if not cfg.include_coverage:
            csv_columns = [col for col in csv_columns if
                           col.rstrip('0123456789') not in cov_columns]

        writer = csv.DictWriter(out_file, fieldnames=csv_columns,
                                extrasaction='ignore')
        writer.writeheader()

        for gene in genes.values():
            gene_dict = gene.__dict__
            to_delete = []
            for attr, value in list(gene_dict.items()):
                if not cfg.include_coverage and attr + '_' in cov_columns:
                    to_delete.append(attr)
                    continue
                if type(value) == dict:
                    if attr == "coords":
                        name = "PC_"
                    else:
                        name = attr + "_"
                    for index, sub_val in value.items():
                        gene_dict[name + str(index)] = str2float_or_int(sub_val,
                                                                        2)
                    to_delete.append(attr)
                elif attr == 'bh_evalue' or attr == 'c_pct_assemby_len':
                    # don't round e-value and assembly length percentage (ratio)
                    pass
                else:
                    gene_dict[attr] = str2float_or_int(value, 2)

            for key in to_delete:
                del gene_dict[key]

            writer.writerow(gene_dict)


def write_tmp_query_info(output_path, genes, queryID, taxon_name):
    """
    Write temporary file which has information about query taxon.

    Args:
      output_path:
      genes:
      queryID:

    Returns:

    """
    query_label = ident_query_label(genes, queryID)
    with open(output_path + 'tmp/tmp.query_label', 'w') as tmp_file:
        tmp_file.write(query_label + '\n')
        tmp_file.write(taxon_name + '\n')


def taxon_sunburst(cfg, genes):
    """ create input for hierarchic taxon sunburst plot

    :param cfg:
    :param genes:
    :return:
    """

    sunburst_dict = {}
    childs_dict = {}  # key: taxon_parent , value: [taxon_childs]
    parent_dict = {}  # key: taxon_child , value: taxon_parent
    leaves_count = {}  # key: leaf taxon, value = count

    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            continue
        g_taxon = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB)
        ta_lineage = g_taxon.name_lineage
        for i in range(len(ta_lineage) - 1):
            if i == 0:
                if ta_lineage[i] in leaves_count.keys():
                    leaves_count[ta_lineage[i]] += 1
                else:
                    leaves_count[ta_lineage[i]] = 1
            if ta_lineage[i + 1] in childs_dict.keys():  # parent
                childs_dict[ta_lineage[i + 1]].add(ta_lineage[i])
            else:
                childs_dict[ta_lineage[i + 1]] = set([ta_lineage[i]])
            parent_dict[ta_lineage[i]] = ta_lineage[i + 1]

    while leaves_count:
        del_leaves = []
        for child, count in leaves_count.items():
            if childs_dict.get(child) != None:
                continue
            del_leaves.append(child)
            sunburst_dict[child] = {"name": child, "value": count}
            parent = parent_dict.get(child)
            while parent != None and child in childs_dict.get(parent):
                if parent in sunburst_dict.keys():
                    sunburst_dict[parent]["children"].append(
                        sunburst_dict[child])
                else:
                    if parent in leaves_count.keys():
                        sunburst_dict[parent] = {"name": parent, "children": [
                                                                                 sunburst_dict[
                                                                                     child]] +
                                                                             [{
                                                                                  "name": '',
                                                                                  "value":
                                                                                      leaves_count[
                                                                                          parent]}]}
                        del_leaves.append(parent)
                    else:
                        sunburst_dict[parent] = {"name": parent, "children": [
                            sunburst_dict[child]]}
                childs_dict[parent].remove(child)
                if len(childs_dict.get(parent)) > 1:
                    break
                child = parent
                parent = parent_dict.get(child)

        for key in del_leaves:
            del leaves_count[key]

    sunburst_dict = sunburst_dict.get("root")

    check_sunburst(genes, sunburst_dict)

    print(json.dumps(sunburst_dict, indent=""))
    return json.dumps(sunburst_dict)


def check_sunburst(genes, dict):
    """ checks the output of taxon_sunburst() """

    leaves_count = {}

    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            continue
        if gene.taxon_assignmentID in leaves_count.keys():
            leaves_count[gene.taxon_assignmentID] += 1
        else:
            leaves_count[gene.taxon_assignmentID] = 1

    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            continue
        g_taxon = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB)
        ta_lineage = g_taxon.name_lineage[::-1]
        tmp_dict = dict
        for taxon in ta_lineage[1:]:
            for child in tmp_dict.get("children"):
                if taxon == child.get("name"):
                    tmp_dict = child
                    bool = True
                    break
                else:
                    bool = False
            if not bool:
                sys.exit("lineage in sunburst dictionary not correct")
        if "children" in tmp_dict.keys():
            for child in tmp_dict["children"]:
                if child["name"] == "":
                    tmp_dict = child
        if tmp_dict.get("value") != leaves_count.get(gene.taxon_assignmentID):
            sys.exit("taxon count in sunburst dictionary not correct")


###############################################################################
###############################################################################


def run_assignment(cfg, pca_coordinates, TAX_DB_local):
    """

    Args:
      cfg:

    Returns:

    """

    global TAX_DB
    TAX_DB = TAX_DB_local

    global missing_taxids
    missing_taxids = set()

    if cfg.assignment_mode == "quick":
        tax_assignment_path_1, tax_assignment_path_2 = cfg.tax_assignment_path
    else:
        tax_assignment_path = cfg.tax_assignment_path
    num_threads = f' -p {cfg.threads} ' \
                    if cfg.threads != 'auto' else ' '

    tmp_prot_path = cfg.output_path + "tmp/tmp.subset.protein.fasta"

    diamond_cmd = f'{cfg.diamond} blastp {num_threads} ' \
                  f'-f 6 qseqid sseqid pident length mismatch gapopen qstart' \
                  f' qend sstart send evalue bitscore staxids sscinames ' \
                  f'-b2.0 --tmpdir /dev/shm --sensitive -c1 ' \
                  f'--top 10 -q "{tmp_prot_path}" -d "{cfg.database_path}"'

    if cfg.taxon_exclude:
        taxon_exclude = f' --taxon-exclude ' \
                        f'"{get_id_for_rank_of_species(cfg.taxon_id, cfg.exclusion_rank)}"'
    else:
        taxon_exclude = ''
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

    # read file
    genes, header = read_raw_genes(cfg)
    header, omitted_genes = assign_pca_coords(genes, pca_coordinates, header)

    prots = prot_gene_matching(cfg.output_path, cfg.gff_path, genes, cfg,
                               omitted_genes)

    perform_diamond = cfg.compute_tax_assignment and not cfg.update_plots

    if perform_diamond:
        subset_prots_longest_cds(genes, cfg.proteins_path, tmp_prot_path)
    if cfg.assignment_mode == 'quick':
        perform_quick_search_1(perform_diamond, diamond_cmd,
                               tax_assignment_path_1,
                               cfg.quick_mode_search_rank,
                               taxon_exclude, cfg.taxon_id)
        taxonomic_assignment(tax_assignment_path_1, genes, prots, cfg.taxon_id)
        perform_quick_search_2(perform_diamond, diamond_cmd,
                               tax_assignment_path_2, cfg.quick_mode_match_rank,
                               genes, tmp_prot_path, taxon_exclude,
                               cfg.taxon_id)
        taxonomic_assignment(tax_assignment_path_2, genes, prots, cfg.taxon_id)
    elif cfg.assignment_mode == 'exhaustive':
        perform_exhaustive_search(perform_diamond, diamond_cmd,
                                  tax_assignment_path, taxon_exclude)
        taxonomic_assignment(tax_assignment_path, genes, prots, cfg.taxon_id)
    else:
        logging.error('Assignment mode not one of quick or exhaustive')
        sys.exit()

    merge_labels(genes, cfg.taxon_id, cfg.merging_labels)
    if num_groups_plot:
        unlabeled_genes, num_labels = filter_labeled_genes(genes)
        iterative_merging_top(unlabeled_genes, num_groups_plot - num_labels)
    else:
        set_unassigned_labels(genes)

    write_tmp_query_info(cfg.output_path, genes, cfg.taxon_id, cfg.taxon_name)

    if missing_taxids:
        logging.info(
            "The following Taxon ID(s) could not be found in the NCBI (skipped in taxonomic assignment):")
        logging.info(missing_taxids)

    pathlib.Path(cfg.output_path + 'taxonomic_assignment/').mkdir(parents=True,
                                                                  exist_ok=True)

    # taxon_sunburst(cfg, genes)

    taxonomy_summary(cfg, genes)
    write_output(cfg, genes, header)
    return genes


def main():
    """ """
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkData.cfg2obj(config_path)

    run_assignment(cfg)


if __name__ == '__main__':
    main()
