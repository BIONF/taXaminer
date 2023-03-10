#!/usr/bin/env python

"""Perform taxonomic assignment

Expects processed config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

from . import checkInput

import taxopy
import sys
import csv
import pathlib
import subprocess
import logging
import pandas as pd


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


def get_id_for_rank_of_species(target_id, rank):
    """
    Returns ID of taxon at given rank for given taxon ID.

    Args:
      target_id(int): NCBI ID of taxon for wich to find ID at given rank
      rank(str): rank of which to get the NCBI ID of

    Returns:
       (int) NCBI ID of taxon at the given rank for given taxon ID
     """
    query_taxon = init_taxopyon(target_id)
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


def init_taxopyon(taxid):
    """ init taxopy taxon object
    """

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
        taxon = init_taxopyon(taxID)
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
        taxon = init_taxopyon(gene.taxon_assignmentID)
        ranks = taxon.rank_name_dictionary
        merger_rank = ranks.get(rank)
        if merger_rank:
            gene.plot_label = merger_rank
            gene.plot_labelID = get_id_for_rank_of_species(
                gene.taxon_assignmentID, rank)
        else:
            gene.plot_label = gene.taxon_assignment
            gene.plot_labelID = gene.taxon_assignmentID


def merge_assignments_at_id(assignments_df, ids):
    """
    Merge taxonomic assignments at given NCBI taxon IDs.

    Args:
      assignments_df(dict): dictionary of genes {g_name: gene object}
      ids(int,list): list or int of taxon ID(s) at which to merge
    """

    if type(ids) != list:
        ids = [ids]

    for g_name, gene in assignments_df.items():
        if gene.taxon_assignmentID == 'NA':
            gene.plot_label = 'Unassigned'
            gene.plot_labelID = 'NA'
            continue

        taxon = init_taxopyon(gene.taxon_assignmentID)
        lineage = taxon.taxid_lineage
        for id in ids:
            if id in lineage:
                id_merger_taxon = init_taxopyon(id)
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
                    lca_candidates.append(init_taxopyon(taxon))
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


def merge_labels(assignments_df, target_id, merging_labels):
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
                query_merger_id = get_id_for_rank_of_species(target_id,
                                                             merger_label)
                break
        else:
            if type(merger_label) == int or merger_label.isdigit():
                merger_ids.append(
                    int(merger_label))  # ncbiID of <merger_label> to merge at

    # if there is an -all merger-label; this overwrites everything
    # else; merge at query_merger_id and merger_ids
    if all_merger_rank:
        merge_assignments_at_rank(assignments_df, all_merger_rank)
    else:
        if query_merger_id:
            merge_assignments_at_id(assignments_df, query_merger_id)
        if merger_ids:
            merge_assignments_at_id(assignments_df, merger_ids)


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


def ident_query_label(assignments_df, target_taxon):
    """
    Identify which label represents the query species.

    Args:
      genes:
      target_id:

    Returns:

    """

    target_lineage = target_taxon.taxid_lineage

    label_candidate = (None, len(target_lineage))
    for assignment in assignments_df['taxon_assignmentID'].unique():
        if assignment in target_lineage:
            if target_lineage.index(assignment) < label_candidate[1]:
                label_candidate = (assignment, target_lineage.index(assignment))

    return label_candidate[0]


###############################################################################
########################## TAXONOMIC ASSIGNMENT ###############################
###############################################################################

def get_subset_id(target_id, subset_marker):
    """
    Return input if input.isdigit() or return ID for rank of query.

    Args:
      target_id:
      subset_marker:

    Returns:

    """
    if type(subset_marker) == int:
        return subset_marker
    elif subset_marker.isdigit():
        return int(subset_marker)
    else:
        return get_id_for_rank_of_species(target_id, subset_marker)


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


def perform_quick_search_1(cfg, perform_diamond, diamond_cmd,
                           quick_mode_search_rank, target_exclude, target_taxon):
    """
    Perform first DIAMOND run for quick assignment mode.

    Args:
      perform_diamond:
      diamond_cmd:
      tax_assignment_path_1:
      quick_mode_search_rank:
      target_exclude:
      target_id:
    """
    logging.info("Quick mode for taxonomic assignment selected")

    if perform_diamond:
        # quick_mode_search_rank can be either taxon ID or rank
        quick_mode_search_id = get_subset_id(target_taxon.taxid, quick_mode_search_rank)

        if target_exclude:
            q1_exclude = target_exclude.rstrip('"') + ',' + ','.join(
                diamond_inclusion_by_exclusion(quick_mode_search_id))
        else:
            q1_exclude = f' --taxon-exclude ' \
                         f'"{",".join(diamond_inclusion_by_exclusion(quick_mode_search_id))}"'
        diamond_cmd_1 = f'{diamond_cmd} -o "{cfg.diamond_results_path}" {q1_exclude}'
        ###taxonlist###
        # diamond_cmd_1 = diamond_cmd + diamond_o1
        # + ' --taxonlist "' + str(quick_mode_search_id)
        # + '" --taxon-exclude "' + str(target_id) + '"'

        logging.debug(f"Diamond command for "
                      f"quick search round 1:\n{diamond_cmd_1}")
        run_diamond(diamond_cmd_1)


def perform_quick_search_2(cfg, perform_diamond, diamond_cmd,
                           quick_mode_match_rank, assignment_df, tmp_prot_path,
                           target_exclude, target_taxon):
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
    """

    if perform_diamond:
        quick_mode_match_id = get_subset_id(target_taxon.taxid, quick_mode_match_rank)

        filter_prots_quick_hits(assignment_df, quick_mode_match_id,
                                tmp_prot_path, tmp_prot_path)

        diamond_cmd_2 = f'{diamond_cmd} -o "{cfg.diamond_results_path}" {target_exclude}'
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
        diamond_cmd = f'{diamond_cmd} -o "{cfg.diamond_results_path}" {target_exclude}'
        logging.debug(
            f"Diamond command for exhaustive search:\n{diamond_cmd}")
        run_diamond(diamond_cmd)


def assess_closest_of_hits(taxIDs, target_taxon):
    """
    for a list of taxa IDs, get the one where LCA with target_id is closest to query

    Args:
      taxIDs:
      target_id:

    Returns:

    """

    target_lineage = target_taxon.taxid_lineage

    # (taxon, index in lineage, lca with target)
    closest_hit = (None, len(target_lineage) + 1)
    for taxID in taxIDs:     
        taxon = init_taxopyon(taxID)
        if taxon:
            lca = taxopy.find_lca([target_taxon, taxon], TAX_DB)
            lca_index = target_lineage.index(lca.taxid)
            # check if new best hit is closer to query in lineage as current best hit
            if closest_hit[1] > lca_index:              
                closest_hit = (taxon, lca_index, lca)

    return closest_hit[0], closest_hit[2]

def compute_lca(taxon_ids):
    """ compute the LCA of a list of taxon ID """

    taxa_list = [init_taxopyon(taxon_id) for taxon_id in taxon_ids]
    if len(taxa_list) > 1:
        return taxopy.find_lca(taxa_list, TAX_DB)
    else:
        return taxa_list[0]


def write_hits2file(cfg, file_path, hits):

    with open(file_path, 'w') as file:
        for hit in hits:
            if cfg.slim_diamond_results:
                file.write('\t'.join(hit[:3]+hit[10:]) + '\n')
            else:
                file.write('\t'.join(hit) + '\n')


def calc_assignment(assignments_df, gene_id, hit_list, target_taxon):
    """
    Calculate the LCA for each gene based on taxonomic hits.

    Args:
      genes:

    Returns:
    :param target_taxon:

    """

    hit_ids = [id[-2].split(';')[0] for id in hit_list]
    lca = compute_lca(hit_ids)

    assignments_df.at[gene_id,'lcaID'] = lca.taxid
    assignments_df.at[gene_id,'lca'] = lca.name

    if lca.taxid in target_taxon.taxid_lineage:
        # refine the LCA with taxon that has clostest LCA with target
        closest_hit, closest_lca = assess_closest_of_hits(hit_ids, target_taxon)
        reference_hit = [hit for hit in hit_list
                         if str(closest_hit.taxid) in hit[-2]][0]
        assignments_df.at[gene_id, 'ref_hitID'] = closest_hit.taxid
        assignments_df.at[gene_id, 'ref_hit'] = closest_hit.name
        assignments_df.at[gene_id, 'ref_evalue'] = reference_hit[10]
        assignments_df.at[gene_id, 'ref_pident'] = reference_hit[2]
        assignments_df.at[gene_id, 'ref_bitscore'] = reference_hit[11]
        assignments_df.at[gene_id, 'refined_lcaID'] = closest_lca.taxid
        assignments_df.at[gene_id, 'refined_lca'] = closest_lca.name
        assignments_df.at[gene_id, 'taxon_assignmentID'] = closest_lca.taxid
        assignments_df.at[gene_id, 'taxon_assignment'] = closest_lca.name

        hit_list.remove(reference_hit)
        return [reference_hit] + hit_list

    else:
        ref_hit_taxon = init_taxopyon(hit_list[0][-2])
        assignments_df.at[gene_id, 'ref_hitID'] = ref_hit_taxon.taxid
        assignments_df.at[gene_id, 'ref_hit'] = ref_hit_taxon.name
        assignments_df.at[gene_id, 'ref_evalue'] = hit_list[0][10]
        assignments_df.at[gene_id, 'ref_pident'] = hit_list[0][2]
        assignments_df.at[gene_id, 'ref_bitscore'] = hit_list[0][11]
        assignments_df.at[gene_id, 'taxon_assignmentID'] = lca.taxid
        assignments_df.at[gene_id, 'taxon_assignment'] = lca.name

        return hit_list



def process_multi_hit(hit, target_taxon):

    split_taxids = hit[-2].split(';')
    lca = compute_lca(split_taxids)

    if lca.taxid in target_taxon.taxid.lineage:
        closest_hit, closest_lca = assess_closest_of_hits(split_taxids, target_taxon)
        return [f"{closest_hit.taxid};LCA:{lca.taxid}",
                f"{closest_hit.name};LCA:{lca.name}"]
    else:
        return [f"{lca.taxid};LCA:{lca.taxid}", f"{lca.name};LCA:{lca.name}"]


def add_hit2list(list, hit, target_taxon):

    if ';' in hit[-2]:  # hit is a multihit
        cropped_hit = process_multi_hit(hit, target_taxon)
        list.append(hit[:-2]+cropped_hit)
    else:
        list.append(hit)


def set_tax_assignments(cfg, tax_assignment_path, assignments_df, target_taxon):
    """
    Read DIAMOND output file and assign hits to genes.

    Args:
      tax_assignment_path:
      prots:
      target_id:

    Returns:

    """

    with open(cfg.diamond_results_path, 'r') as diamond_hits:
        # cols: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore taxid taxname
        # tax assignments list is sorted by query sequence
        # thus, all hits for one gene are following one another

        first_spline = next(diamond_hits).strip().split('\t')
        gene = assignments_df.loc[
            assignments_df['transcript_id'] == first_spline[0]]
        gene_assignments = []
        add_hit2list(gene_assignments, first_spline, target_taxon)
        for line in diamond_hits:
            # first is "best hit" in terms of alignment score
            spline = line.strip().split('\t')
            if spline[0] == gene.transcript_id.item():
                add_hit2list(gene_assignments, spline, target_taxon)
            else:
                # compute the assignment with list of hits
                if not gene.empty:
                    hits = calc_assignment(assignments_df, gene.index.item(),
                                           gene_assignments, target_taxon)
                    write_hits2file(cfg, tax_assignment_path, hits)

                # new gene
                gene = assignments_df.loc[
                    assignments_df['transcript_id'] == spline[0]]
                gene_assignments = []
                add_hit2list(gene_assignments, spline, target_taxon)

        if not gene.empty:
            hits = calc_assignment(assignments_df, gene.index.item(), gene_assignments,
                            target_taxon)
            write_hits2file(cfg, tax_assignment_path, hits)


def taxonomic_assignment(cfg, gff_df, assignments_df, tax_assignment_path, target_taxon):
    """
    Run functions to save taxonomic assignment to gene objects

    Args:
      tax_assignment_path:
      genes:
      prots:
      target_id:
    """
    set_tax_assignments(cfg, tax_assignment_path, assignments_df, target_taxon)

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
    longest_transcripts = gff_df.loc[gff_df['type'] == 'gene', 'transcript_id'].to_list()

    logging.info(
        f"{len(longest_transcripts)} proteins written to fasta file for "
        f"taxonomic assignment (subsetting for longest CDS)")
    subset_protein_fasta(proteins_path, longest_transcripts, path_out,
                         "include")


def detect_appropriate_assignment(gene, quick_mode_match_id):


    tax_lineage = init_taxopyon(gene.taxon_assignmentID).taxid_lineage
    # if quick_mode_match_id not in tax_lineage or quick_mode_match_id == tax_lineage[0]: # <- option 2: exclusive
    if quick_mode_match_id not in tax_lineage:  # <- option 1: inclusive
        # inclusive: if match id equals class, assignment must be at least same class as query
        # exclusive: if match id equals class, assignment must be at least one rank closer to query
        return gene.transcript_id
    else:
        return None




def filter_prots_quick_hits(assignments_df, quick_mode_match_id, prot_path_in,
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

    no_match = assignments_df.apply(
        lambda row: detect_appropriate_assignment(row, quick_mode_match_id),
        axis=1).to_list()

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

    include_taxon = init_taxopyon(include_id)
    include_lineage = include_taxon.taxid_lineage
    for parent in include_lineage[1:]:
        childs = get_children_ids(parent)
        exclude_list = exclude_list + [str(child) for child in childs if
                                       child not in include_lineage]

    return exclude_list


def taxonomy_summary(cfg, assignments_df, target_taxon, query_label):
    """
    computing a score to assess the quality of the assembly

    Args:
      cfg:
      assignments_df:

    Returns:

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
        g_taxon = init_taxopyon(row['taxon_assignmentID'])
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

def assignment_lineages(assignments_df):

    rank_subset = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                   'family', 'genus', 'species', 'strain']

    for assignment in assignments_df['taxon_assignmentID'].unique():
        assigned_taxon = init_taxopyon(assignment)
        ranks_of_interest = {key: assigned_taxon.rank_name_dictionary.get(key)
                             for key in rank_subset}


###############################################################################
###############################################################################


def run_assignment(cfg, gff_df, pca_coordinates, TAX_DB_local):
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


    tmp_prot_path = cfg.output_path + "tmp/tmp.subset.protein.fasta"

    diamond_cmd = f'{cfg.diamond} blastp -p {cfg.threads} ' \
                  f'-f 6 qseqid sseqid pident length mismatch gapopen qstart' \
                  f' qend sstart send evalue bitscore staxids sscinames ' \
                  f'-b2.0 --tmpdir /dev/shm --sensitive -c1 ' \
                  f'--top 10 -q "{tmp_prot_path}" -d "{cfg.database_path}"'

    target_taxon = init_taxopyon(cfg.taxon_id)

    if cfg.target_exclude:
        target_exclude = f' --taxon-exclude ' \
                        f'"{get_id_for_rank_of_species(target_taxon.taxid, cfg.exclusion_rank)}"'
    else:
        target_exclude = ''


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

    perform_diamond = cfg.compute_tax_assignment and not cfg.update_plots

    # init empty dataframe for taxonomic assignments
    num_genes = gff_df.loc[gff_df['type'] == 'gene'].shape[0]
    assignments_dict = {
        'transcript_id': gff_df.loc[gff_df['type'] == 'gene', 'transcript_id'].to_list(),
        'ref_hit': [None for i in range(num_genes)],
        'ref_hitID': [None for i in range(num_genes)],
        'ref_pident': [None for i in range(num_genes)],
        'ref_evalue': [None for i in range(num_genes)],
        'ref_bitscore': [None for i in range(num_genes)],
        'lca': [None for i in range(num_genes)],
        'lcaID': [None for i in range(num_genes)],
        'refined_lca': [None for i in range(num_genes)],
        'refined_lcaID': [None for i in range(num_genes)],
        'taxon_assignment': [None for i in range(num_genes)],
        'taxon_assignmentID': [None for i in range(num_genes)]
    }
    assignments_df = pd.DataFrame(assignments_dict,
        index = gff_df.loc[gff_df['type'] == 'gene'].index)

    if perform_diamond:
        subset_prots_longest_cds(gff_df, cfg.proteins_path, tmp_prot_path)
    if cfg.assignment_mode == 'quick':
        perform_quick_search_1(cfg, perform_diamond, diamond_cmd, cfg.quick_mode_search_rank,
                               target_exclude, cfg.taxon_id)
        taxonomic_assignment(cfg, gff_df, assignments_df,
                             tax_assignment_path_1, target_taxon)
        perform_quick_search_2(cfg, perform_diamond, diamond_cmd,
                               cfg.quick_mode_match_rank,
                               assignments_df, tmp_prot_path, target_exclude,
                               target_taxon)
        taxonomic_assignment(cfg, gff_df, assignments_df,
                             tax_assignment_path_2, target_taxon)
    elif cfg.assignment_mode == 'exhaustive':
        perform_exhaustive_search(cfg, perform_diamond, diamond_cmd, target_exclude)
        taxonomic_assignment(cfg, gff_df, assignments_df,
                             tax_assignment_path, target_taxon)
    else:
        logging.error('Assignment mode not one of quick or exhaustive')
        sys.exit()


    # merge_labels(assignments_df, target_taxon, cfg.merging_labels)
    # if num_groups_plot:
    #     unlabeled_genes, num_labels = filter_labeled_genes(assignments_df)
    #     iterative_merging_top(unlabeled_genes, num_groups_plot - num_labels)
    # else:
    #     set_unassigned_labels(assignments_df)
    #
    query_label = ident_query_label(assignments_df, target_taxon)


    if missing_taxids:
        logging.info(
            "The following Taxon ID(s) could not be found in the NCBI (skipped in taxonomic assignment):")
        logging.info(missing_taxids)

    pathlib.Path(cfg.output_path + 'taxonomic_assignment/').mkdir(parents=True,
                                                                  exist_ok=True)

    assignment_lineages(assignments_df)
    taxonomy_summary(cfg, assignments_df, target_taxon, query_label)
    write_output(cfg, assignments_df)
    return assignments_df, query_label


def main():
    """ """
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    run_assignment(gff_df, cfg)


if __name__ == '__main__':
    main()
