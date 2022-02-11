# -*- coding: utf-8 -*-

# reads FASTA, GFF3 and per_base_coverage.txt file and produces gene_info directory
# author: Simonida Zehr and Freya Arthen
# date: 18 June 2019

import numpy as np
import pathlib # to create directories
import operator # for quick comparisons
import scipy.stats as stats # for Pearson's R
from itertools import product as itertools_product # to generate all possible oligonucleotides from base alphabet
from Bio.Seq import Seq as BioPython_Seq # to count oligonucleotides (also overlapping ones! not implemented in normal count)
import sys # parse command line arguments
import logging

from classes import Assembly, Contig, Gene
import prepare_and_check

# ====================== CONSTANTS ======================
# define nucleotide alphabet
NUC_ALPHABET = ['A', 'C', 'G', 'T']

# generate an array holding all possible (4^4 = 256) tetranucleotides as strings
TETRANUCLEOTIDES = [''.join(i) for i in itertools_product(NUC_ALPHABET, repeat = 4)]
# generate an array holding all possible (4^3 = 64) trinucleotides as strings
TRINUCLEOTIDES = [''.join(i) for i in itertools_product(NUC_ALPHABET, repeat = 3)]
# generate an array holding all possible (4^2 = 16) dinucleotides as strings
DINUCLEOTIDES = [''.join(i) for i in itertools_product(NUC_ALPHABET, repeat = 2)]

# ==============================================================
# ========================= FUNCTIONS ==========================
# ==============================================================
def compare(given, relation, crit):
    """
    Returns True, if the given value relates to the critical value
    in the way that is described by relation.
    i.e. if relation=">" and crit=3, then
    True is returned e.g. if given=2, False is returned if given=4.

    Args:
      given:
      relation:
      crit:

    Returns:

    """
    operators = {">" : operator.gt, "above" : operator.gt,
                 "<" : operator.lt, "below" : operator.lt,
                ">=" : operator.ge, "<=" : operator.le,
                "==" : operator.eq, "eq" : operator.eq}
    return operators[relation](given, crit)


def round_if_float(x, y):
    """
    Returns rounded x (if it is a float).

    Args:
      x: value to be rounded
      y: digits after the decimal point

    Returns:

    """

    if isinstance(x, float) and not np.isnan(x):
        return round(x,y)
    else:
        return x


def invert(x, min_value, max_value):
    """
    Invert the given value on its respective scale,
    as defined by the minimum and the maximum value of the scale
    (if x is not NaN).

    Args:
      x:
      min_value:
      max_value:

    Returns:

    """

    if np.isnan(x):
        return x
    else:
        return ((max_value + min_value) - x)


def rescale(x, old_min, old_max, new_min, new_max, nansub):
    """
    Assigns the given value (x), which is located on a scale
    ranging from old_min to old_max, to a corresponding new value y,
    which is located on a scale ranging from new_min to new_max.
    NaN values are imputed with the value given in the variable nansub.

    Note: The current implementation features a normalisation to range [0,1]
    If you wish to change this, instructions are given below.

    Args:
      x:
      old_min:
      old_max:
      new_min:
      new_max:
      nansub:

    Returns:

    """


    if np.isnan(x):
        y = ((nansub - old_min)/(old_max - old_min))
    else:
        y = ((x - old_min)/(old_max - old_min))
    return y


def remove_prefix(text, prefix):
    """

    Args:
      text:
      prefix:

    Returns:

    """
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def strip_ID(id):
    """
    Remove GFFs prefixes from IDs

    Args:
      id:

    Returns:

    """
    for prefix in ['gene:','gene-','transcript:','transcript-','rna:','rna-','cds:','cds-']:
        id = remove_prefix(id,prefix)
    return id


def impute_array(a, data_array, include_coverage):
    """
    This function takes a data_array, i.e. an array of arrays describing a
    matrix where each row represents a gene and each column a variable.
    NaN values in a variable are imputed with the mean of this variable
    and all values are rescaled to the range [0,1].

    Args:
      a:
      data_array:
      include_coverage:

    Returns:

    """

    # define the indices of the variables to be transformed
    index_of = {
        "c_genelensd" : 6,
        "c_genecovsd" : 11, # e.g. 12th column in the matrix holds c_genecovsd
        "g_lendev_c" : 17,
        "g_covdev_c" : 24,
        "g_gcdev_c" : 31
    }

    # we're going to loop over all genes and store all of the
    # variables to be imputed
    # --> initialise these arrays
    c_genecovsds = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()}
    c_genelensds = []
    g_covdev_cs = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()}
    g_lendev_cs = []
    g_gcdev_cs = []

    # for each gene array in the given matrix:
    for gene_array in data_array:
        for pbc_index in a.get_pbc_paths().keys():
            # get all gene coverage SDs which are not NaNs
            if not np.isnan(gene_array[index_of["c_genecovsd"]][pbc_index]):          # (expl.: index_of["c_genecovsd"] returns 11
                # and add them to the respective array                              # gene_array[11] will give the c_genecovsd value
                c_genecovsds[pbc_index].append(gene_array[index_of["c_genecovsd"]][pbc_index]) # for the current gene)
            # get all gene deviations from contig mean cov
            if not np.isnan(gene_array[index_of["g_covdev_c"]][pbc_index]):
                g_covdev_cs[pbc_index].append(gene_array[index_of["g_covdev_c"]][pbc_index])

        # get all gene length SDs
        if not np.isnan(gene_array[index_of["c_genelensd"]]):
            c_genelensds.append(gene_array[index_of["c_genelensd"]])
        # get all gene deviations from contig mean GCss
        if not np.isnan(gene_array[index_of["g_gcdev_c"]]):
            g_gcdev_cs.append(gene_array[index_of["g_gcdev_c"]])
        # get all gene deviations from contig mean lengths
        if not np.isnan(gene_array[index_of["g_lendev_c"]]):
            g_lendev_cs.append(gene_array[index_of["g_lendev_c"]])
    # --> at the end of this for loop, we have populated the 5 arrays with the (not-NaN) values of the genes
    # thus, each array represents all (not-NaN) values in one variable

    # for each "variable"
    # get minimum and maximum of all values (needed for inversion and rescaling (see below))
    min_c_genecovsds = {pbc_index: min(coverages, default=np.nan) for pbc_index, coverages in c_genecovsds.items()}
    max_c_genecovsds = {pbc_index: max(coverages, default=np.nan) for pbc_index, coverages in c_genecovsds.items()}

    min_c_genelensds = min(c_genelensds, default=np.nan)
    max_c_genelensds = max(c_genelensds, default=np.nan)

    min_g_covdev_cs = {pbc_index: min(coverages, default=np.nan) for pbc_index, coverages in g_covdev_cs.items()}
    max_g_covdev_cs = {pbc_index: max(coverages, default=np.nan) for pbc_index, coverages in g_covdev_cs.items()}

    min_g_lendev_cs = min(g_lendev_cs, default=np.nan)
    max_g_lendev_cs = max(g_lendev_cs, default=np.nan)

    min_g_gcdev_cs = min(g_gcdev_cs, default=np.nan)
    max_g_gcdev_cs = max(g_gcdev_cs, default=np.nan)

    # get the respective mean of each variable (ignoring NaNs and checking for empty lists)
    my_mean1 = {pbc_index: (np.nanmean(coverages) if coverages else np.nan) for pbc_index, coverages in c_genecovsds.items()}
    my_mean2 = np.nanmean(c_genelensds) if g_lendev_cs else np.nan
    my_mean3 = {pbc_index: (np.nanmean(coverages) if coverages else np.nan) for pbc_index, coverages in g_covdev_cs.items()}
    my_mean4 = np.nanmean(g_gcdev_cs) if g_lendev_cs else np.nan
    my_mean5 = np.nanmean(g_lendev_cs) if g_lendev_cs else np.nan

    # loop over all genes, impute NaNs with mean of variable and rescale to given range (here: [0,1])
    for gene_array in data_array:
        if include_coverage:
            for pbc_index in a.get_pbc_paths().keys():
                gene_array[index_of["c_genecovsd"]][pbc_index] = rescale(gene_array[index_of["c_genecovsd"]][pbc_index], \
                                min_c_genecovsds[pbc_index], max_c_genecovsds[pbc_index], 0, 1, my_mean1[pbc_index]) #1.0, 10.0
                gene_array[index_of["g_covdev_c"]][pbc_index] = rescale(gene_array[index_of["g_covdev_c"]][pbc_index], \
                                min_g_covdev_cs[pbc_index], max_g_covdev_cs[pbc_index], 0, 1, my_mean3[pbc_index]) #1.0, 10.0

        gene_array[index_of["c_genelensd"]] = \
        rescale(gene_array[index_of["c_genelensd"]], min_c_genelensds, max_c_genelensds, 0, 1, my_mean2) #1.0, 10.0
        gene_array[index_of["g_gcdev_c"]] = \
        rescale(gene_array[index_of["g_gcdev_c"]], min_g_gcdev_cs, max_g_gcdev_cs, 0, 1, my_mean4) #1.0, 10.0
        gene_array[index_of["g_lendev_c"]] = \
        rescale(gene_array[index_of["g_lendev_c"]], min_g_lendev_cs, max_g_lendev_cs, 0, 1, my_mean5) #1.0, 10.0

    return data_array


def get_lgth_corr_stats(values, lengths):
    """
    Returns mean and standard deviation of the given list of values
    (e.g. coverages, GC contents, ...),
    weighted with the given list of corresponding lengths.
    * values should be a list holding the values (cov, GC, ...) of the sequences
    of interest (e.g. all genes of one contig, or all contigs in the set)
    * lengths should be a list holding the sequence lengths
    The corresponding elements should have the same indices.

    Args:
      values:
      lengths:

    Returns:

    """

    # if weights sum up to 0, weighted average can not be computed
    if sum(lengths) == 0:
        return np.nan, np.nan
    elif all(x is None for x in values):
        return np.nan, np.nan

    # compute weighted mean of values
    weighted_mean = np.average(values, weights = lengths)

    # compute weighted variance
    sum_of_squared_deviations = 0
    for (value, length) in zip(values, lengths):
        sum_of_squared_deviations += length*((value - weighted_mean)**2)
    weighted_variance = sum_of_squared_deviations / np.sum(lengths)

    # compute weighted standard deviation
    weighted_sd = weighted_variance**0.5

    return weighted_mean, weighted_sd


def get_lgth_corr_stats_weighted_by_sq_lgth(values, lengths):
    """
    Returns mean and standard deviation of the given list of values
    (e.g. coverages, GC contents, ...),
    weighted by the squared corresponding lengths of the sequences.
    * values should be a list holding the values (cov, GC, ...) of the sequences
    of interest (e.g. all genes of one contig, or all contigs in the set)
    * lengths should be a list holding the sequence lengths
    The corresponding elements should have the same indices.

    Args:
      values:
      lengths:

    Returns:

    """

    # compute weighted mean of coverages
    squared_lengths = [length**2 for length in lengths]
    weighted_mean = np.average(values, weights = squared_lengths)

    # compute weighted variance
    sum_of_squared_deviations = 0
    for (value, length) in zip(values, squared_lengths):
        sum_of_squared_deviations += length*((value - weighted_mean)**2)
    weighted_variance = sum_of_squared_deviations / np.sum(squared_lengths)

    # compute weighted standard deviation
    weighted_sd = weighted_variance**0.5

    return weighted_mean, weighted_sd


def get_sum_of_contig_lengths(a):
    """
    Returns the sum of the lengths of contigs in contig_list,
    e.g. if these are all contigs originally contained in the assembly,
    then this means that the total assembly length is returned).

    Args:
      a:

    Returns:

    """

    contig_lengths = []
    for contig in a.get_contigs().values():
        contig_lengths.append(contig.get_length())

    return sum(contig_lengths)


def get_Nx(a, x, contig_lengths):
    """
    Returns Nx of given contig set, e.g. if x=50 the N50.

    Args:
      a:
      x:
      contig_lengths:

    Returns:

    """

    # either contig_lengths is a specific list given
    # or no list is given and all contigs are input,
    # then generate list:
    if contig_lengths == "all":
        contig_lengths = []
        for contig in a.get_contigs().values():
            tmp_contig_length = contig.get_length()
            contig_lengths.append(tmp_contig_length)

    # if contig_lengths is a given list,
    # create a new list from it as to not manipulate the original one
    my_contig_lengths = list(contig_lengths)

    # --------- Calculate Nx ---------
    # sort in descending order
    my_contig_lengths.sort(reverse=True)
    # get total length of contig set
    total_length = sum(my_contig_lengths)
    # determine which percentage
    percentage = x / 100
    # determine minimum length contig is required to have
    min_length = total_length * percentage

    # start summing contig lengths until X% of total assembly length is reached
    cumulative_sum = 0
    for length in my_contig_lengths:
        cumulative_sum += length
        if cumulative_sum >= min_length:
            return length


def compute_coverage_and_positional_info(a):
    """Uses per base coverage information to compute & set
    coverage mean and SD for all contigs and genes.
    Computes & sets info on gene length mean & SD.
    Computes & sets absolute positions of genes.
    Computes & sets info on gene GC mean & SD.
    Set info about single exon genes.

    Args:
      a:

    Returns:

    """

    all_genes_coverage = []

    # for ALL contigs (with coverage info)
    for contig in a.get_contigs().values():

        # compute coverage info for this contig:
        contig.compute_own_coverage_info()

        #------Length, GC & Coverage Info Genes-------
        # if contig is geneless
        # set its gene coverage mean & SD to np.nan
        if contig.geneless_info() == 1:
            contig.set_gene_coverage_mean_and_sd(
                            {pbc_index: np.nan for pbc_index in a.get_pbc_paths().keys()},
                            {pbc_index: np.nan for pbc_index in a.get_pbc_paths().keys()})
            contig.set_gene_lengths_mean_and_sd(np.nan, np.nan)

        else:
            # compute coverage info for the genes on this contig
            #  (also set length_of_covered_bases for each gene)
            contig.compute_gene_coverage_info(a)
            # compute info on absolute positions of genes
            contig.compute_gene_positional_info(a)

            # to store all ACTUAL gene lengths on this contig
            contig_gene_lengths = []
            # to store how many bases are covered in each gene on this contig
            contig_gene_covered_base_lengths = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()}
            # to store each gene's mean GC value
            contig_gene_gcs = []
            # to store each gene's mean coverage
            contig_gene_covs = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()}

            for gene_name in contig.get_genes():
                gene = a.get_gene(gene_name)

                contig_gene_lengths.append(gene.get_length())
                contig_gene_gcs.append(gene.get_gc_content())
                for pbc_index in a.get_pbc_paths().keys():
                    contig_gene_covered_base_lengths[pbc_index].append(gene.get_length_of_covered_bases(pbc_index))
                    contig_gene_covs[pbc_index].append(gene.get_coverage(pbc_index))


            # lengths do not have to be length corrected

            # compute length corrected gene coverage mean and SD for this contig
            contig_weighted_gene_cov_mean = {}
            contig_weighted_gene_cov_sd = {}
            for pbc_index in a.get_pbc_paths().keys(): #for each coverage file
                cov_mean, cov_sd = get_lgth_corr_stats( \
                    contig_gene_covs[pbc_index], \
                    contig_gene_covered_base_lengths[pbc_index])
                contig_weighted_gene_cov_mean[pbc_index] = cov_mean
                contig_weighted_gene_cov_sd[pbc_index] = cov_sd

            # compute length corrected gene GC mean and SD for this contig
            contig_weighted_gene_gc_mean, contig_weighted_gene_gc_sd \
            = get_lgth_corr_stats(contig_gene_gcs, contig_gene_lengths)

            # if there is only one gene on this contig
            if len(contig.get_genes()) == 1:
                # set length corrected mean as computed, but SD to NaN
                contig.set_gene_lengths_mean_and_sd(np.mean(contig_gene_lengths),
                                np.nan) # set standard deviation to NaN
                contig.set_gene_coverage_mean_and_sd(contig_weighted_gene_cov_mean,
                                {pbc_index: np.nan for pbc_index in a.get_pbc_paths().keys()})
                                # set standard deviation to NaN
                contig.set_gene_gc_mean_and_sd(contig_weighted_gene_gc_mean,
                                              np.nan)


            # if there is more than one gene, also set standard deviation
            else:
                # compute & set info on gene length mean and SD in contig
                contig.set_gene_lengths_mean_and_sd(np.mean(contig_gene_lengths),
                                                np.std(contig_gene_lengths))

                # set info on weighted gene cov mean and SD on contig:
                contig.set_gene_coverage_mean_and_sd(contig_weighted_gene_cov_mean,
                                                     contig_weighted_gene_cov_sd)

                # set info on weighted gene GC mean and SD on contig:
                contig.set_gene_gc_mean_and_sd(contig_weighted_gene_gc_mean,
                                              contig_weighted_gene_gc_sd)


# -------- START: FUNCTIONS FOR TETRANUCLEOTIDE Z-SCORE / CORRELATION COEFF -----


def get_observed_oligonuc_freq(given_seq, oligonuc):
    """
    Given a sequence and a oligonucleotide of interest

    Args:
      given_seq:
      oligonuc:

    Returns:

    """
    return BioPython_Seq(given_seq).count_overlap(oligonuc)

#TODO: merge the next three functions into one
def get_tetranuc_freqs(given_seq):
    """
    Returns dictionary mapping each of the 4^4 = 256 possible tetranucleotides
    to its observed frequency in the given sequence.

    Args:
      given_seq:

    Returns:

    """
    return {tetranuc : get_observed_oligonuc_freq(given_seq, tetranuc) for tetranuc in TETRANUCLEOTIDES}

def get_trinuc_freqs(given_seq):
    """
    Returns a dictionary mapping each of the 4^3 = 64 possible trinucleotides
    to its observed frequency in the given sequence.

    Args:
      given_seq:

    Returns:

    """
    return {trinuc : get_observed_oligonuc_freq(given_seq, trinuc) for trinuc in TRINUCLEOTIDES}


def get_dinuc_freqs(given_seq):
    """
    Returns a dictionary mapping each of the 4^2 = 16 possible dinucleotides
    to its observed frequency in the given sequence.

    Args:
      given_seq:

    Returns:

    """
    return {dinuc : get_observed_oligonuc_freq(given_seq, dinuc) for dinuc in DINUCLEOTIDES}


def update_observed_oligofreqs(old_observed_tetranuc_freqs, old_observed_trinuc_freqs, old_observed_dinuc_freqs,
                            additional_tetranuc_freqs, additional_trinuc_freqs, additional_dinuc_freqs):
    """
    Takes six input dictionaries (string to integer mapping),
    the first 3 mapping all possible tetranucleotides/trinucleotides/dinucleotides
    to their respective frequencies observed in the dataset parsed SO FAR (total count)
    and last 3 mapping all possible tetranucleotides/trinucleotides/dinucleotides
    to their respective frequencies observed in the added sequences (additional count).

    Returns 3 dictionaries, mapping the tetra-/tri-/dinucleotides
    to their updated frequencies (new total counts).

    Args:
      old_observed_tetranuc_freqs:
      old_observed_trinuc_freqs:
      old_observed_dinuc_freqs:
      additional_tetranuc_freqs:
      additional_trinuc_freqs:
      additional_dinuc_freqs:

    Returns:

    """

    updated_observed_tetranucs = {tetranuc : (old_observed_tetranuc_freqs[tetranuc]
                                + additional_tetranuc_freqs[tetranuc])
                                for tetranuc in TETRANUCLEOTIDES}
    updated_observed_trinucs = {trinuc : (old_observed_trinuc_freqs[trinuc]
                                + additional_trinuc_freqs[trinuc])
                                for trinuc in TRINUCLEOTIDES}
    updated_observed_dinucs = {dinuc : (old_observed_dinuc_freqs[dinuc]
                                + additional_dinuc_freqs[dinuc])
                                for dinuc in DINUCLEOTIDES}

    return updated_observed_tetranucs, updated_observed_trinucs, updated_observed_dinucs


def get_expected_tetranuc_freq(tetranuc, trinuc_freqs, dinuc_freqs):
    """
    Returns the expected frequency for the given tetranucleotide.

    Args:
      tetranuc:
      trinuc_freqs:
      dinuc_freqs:

    Returns:

    """
    # get prefix of len 3 of tetranucleotide (n1n2n3)
    first3bases = tetranuc[:3]
    # and corresponding frequency
    obsfreq_first3bases = trinuc_freqs[first3bases]
    # get suffix of len 3 of tetranucleotide (n2n3n4)
    last3bases = tetranuc[-3:]
    # and corresponding frequency
    obsfreq_last3bases = trinuc_freqs[last3bases]

    # get the two bases in the middle of the tetranuc (n2n3)
    middle2bases = tetranuc[1:-1]
    # and corresponding frequency
    obsfreq_middle2bases = dinuc_freqs[middle2bases]

    # to prevent division by 0:
    if obsfreq_middle2bases == 0:
        return np.nan
    else:
        expected_tetranuc_freq = \
        (obsfreq_first3bases*obsfreq_last3bases)/obsfreq_middle2bases
        return expected_tetranuc_freq


def get_approximate_variance_tetranuc(tetranuc, expfreq_tetranuc, trinuc_freqs, dinuc_freqs):
    """
    Returns an approximation of the variance of the
    observed frequency of the given tetranucleotide.

    Args:
      tetranuc:
      expfreq_tetranuc:
      trinuc_freqs:
      dinuc_freqs:

    Returns:

    """

    # get prefix of len 3 of tetranucleotide (n1n2n3)
    first3bases = tetranuc[:3]
    # and corresponding frequency
    obsfreq_first3bases = trinuc_freqs[first3bases]
    # get suffix of len 3 of tetranucleotide (n2n3n4)
    last3bases = tetranuc[-3:]
    # and corresponding frequency
    obsfreq_last3bases = trinuc_freqs[last3bases]
    # get the two bases in the middle of the tetranuc (n2n3)
    middle2bases = tetranuc[1:-1]
    # and corresponding frequency
    obsfreq_middle2bases = dinuc_freqs[middle2bases]

    # to prevent division by 0:
    if obsfreq_middle2bases == 0:
        return np.nan
    else:
        approx_variance = (expfreq_tetranuc *
                    (((obsfreq_middle2bases - obsfreq_first3bases)
                      *(obsfreq_middle2bases - obsfreq_last3bases))
                     /(obsfreq_middle2bases**2)))

        return approx_variance


def get_z_score(obsfreq_tetranuc, expfreq_tetranuc, approxvar_tetranuc):
    """
    Given the observed frequency of a given tetranucleotide in a given sequence,
    the corresponding expected frequency and its variance,
    this function returns a z-score for this specific tetranucleotide
    in the given sequence.

    Args:
      obsfreq_tetranuc:
      expfreq_tetranuc:
      approxvar_tetranuc:

    Returns:

    """

    if np.isnan(expfreq_tetranuc) or np.isnan(approxvar_tetranuc) or approxvar_tetranuc == 0:
        z_score = np.nan
    else:
        z_score = (obsfreq_tetranuc - expfreq_tetranuc)/(approxvar_tetranuc**0.5)

    return z_score


def get_all_obsfreqs(given_seq):
    """
    Returns dictionaries mapping
    - all possible 256 tetranucleotides
    - all possible 64 trinucleotides
    - all possible 16 dinucleotides
     to their frequencies observed in the given sequence
     (string to integer mapping).

    Args:
      given_seq:

    Returns:

    """
    tetranuc_freqs = get_tetranuc_freqs(given_seq)
    trinuc_freqs = get_trinuc_freqs(given_seq)
    dinuc_freqs = get_dinuc_freqs(given_seq)

    return tetranuc_freqs, trinuc_freqs, dinuc_freqs


def calculate_z_score_vector(tetranuc_freqs, trinuc_freqs, dinuc_freqs):
    """
    Given all tetranucleotide-to-frequency, trinucleotide-to-frequency
    and dinucleotide-to-frequency mappings (dictionaries with key:string, value:int)
    observed in a given sequence, this function returns an array that represents
    the z-scores of all possible 256 tetranucleotides for this sequence.

    Args:
      tetranuc_freqs:
      trinuc_freqs:
      dinuc_freqs:

    Returns:

    """
    z_scores = []
    for tetranuc in TETRANUCLEOTIDES:

        obsfreq_tetranuc = tetranuc_freqs[tetranuc]
        expfreq_tetranuc = get_expected_tetranuc_freq(tetranuc, trinuc_freqs, dinuc_freqs)
        approxvar_tetranuc = get_approximate_variance_tetranuc(tetranuc, expfreq_tetranuc,
                                                            trinuc_freqs, dinuc_freqs)

        z_score = get_z_score(obsfreq_tetranuc, expfreq_tetranuc, approxvar_tetranuc)
        z_scores.append(z_score)

    return z_scores


def pearsonr_with_nans_omitted(vector1, vector2):
    """
    Whenever at least one of the given z-score vectors contains a NaN at an index,
    this element is removed, as well as the corresponding element in the other vector.
    Correlation of these NaN-free vectors is computed. Returns array holding Pearson's R and p-value.

    Args:
      vector1:
      vector2:

    Returns:

    """

    nanfree_vector1 = []
    nanfree_vector2 = []

    for (i, j) in zip(vector1, vector2):
        if not (np.isnan(i) or np.isnan(j)):
            nanfree_vector1.append(i)
            nanfree_vector2.append(j)

    # pearsonr can only be computed if vectors are at least 3 elements long
    # nanfree_vectors are too short when gene has repetive sequence,
    # because the z-score vector is full of nans then
    if len(nanfree_vector1) < 2:
        return [np.nan,np.nan]
    else:
        return stats.pearsonr(nanfree_vector1, nanfree_vector2)


# -------- END: FUNCTIONS FOR TETRANUCLEOTIDE Z-SCORE / CORRELATION COEFF -------


def compute_stats(a, lgth_corr_function):
    """
    Computes and returns basic positional, length, coverage and GC info.
    This function MUST be run in order to obtain cov, GC & length statistics!
    (i.e. mean and SD over all contigs and genes)

    Args:
      a:
      lgth_corr_function:

    Returns:

    """

    contig_gcs = [] # gather mean GC value for each contig
    contig_coverages = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()} # gather mean coverage for each contig
    contig_lengths = [] # gather all ACTUAL contig lengths
    # gather how many bases are COVERED for each contig
    contig_covered_bases_lengths = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()}

    # in the considered assembly part, initialise all possible
    # tetranucleotides, trinucleotides, dinucleotides with frequency 0
    considered_tetranuc_freqs = {nuc : 0 for nuc in TETRANUCLEOTIDES}
    considered_trinuc_freqs = {nuc : 0 for nuc in TRINUCLEOTIDES}
    considered_dinuc_freqs = {nuc : 0 for nuc in DINUCLEOTIDES}

    gene_gcs = [] # gather mean GC value for each gene
    gene_coverages = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()} # gather mean coverage for each gene
    gene_lengths = [] # gather all ACTUAL gene lengths
    # gather how many bases are COVERED for each gene
    gene_covered_bases_lengths = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()}

    for contig_name, contig in a.get_contigs().items():

        if contig_name in a.get_geneless_contigs().keys():
            continue

        #------ Gather Length, Coverage & GC Info Contig ----------------
        # get contig length + add it to list of all contig lengths
        contig_length = contig.get_length()
        contig_lengths.append(contig_length)
        for pbc_index in a.get_pbc_paths().keys():
            # get number of covered bases for this contig + add it to the list
            contig_covered_bases_length = contig.get_length_of_covered_bases(pbc_index)
            contig_covered_bases_lengths[pbc_index].append(contig_covered_bases_length)
            # get contig cov + add it to list of all contig covs
            contig_cov = contig.get_coverage(pbc_index)
            contig_coverages[pbc_index].append(contig_cov)
        # get contig GC + add it to list of all contig GC contents
        contig_gc = contig.get_gc_content()
        contig_gcs.append(contig_gc)

        # for this contig,  initialise all possible
        # tetranucleotides, trinucleotides, dinucleotides with frequency 0
        contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs = \
        contig.get_oligofreqs()

        # update the (considered) overall frequencies with the info from this contig
        considered_tetranuc_freqs, considered_trinuc_freqs, considered_dinuc_freqs = \
        update_observed_oligofreqs(considered_tetranuc_freqs, considered_trinuc_freqs, considered_dinuc_freqs,
                                  contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)

        # ------ Gather Length, Coverage & GC Info Genes -----------------------
        for gene_name in contig.get_genes():
            gene = a.get_gene(gene_name)

            # add lengths of genes on this contig to all gene lengths
            gene_length = gene.get_length()
            gene_lengths.append(gene_length)
            for pbc_index in a.get_pbc_paths().keys():
                # get number of covered bases for this gene + add it to the list
                gene_covered_bases_length = gene.get_length_of_covered_bases(pbc_index)
                gene_covered_bases_lengths[pbc_index].append(gene_covered_bases_length)
                # add covs of genes on this contig to all gene covs
                gene_cov = gene.get_coverage(pbc_index)
                gene_coverages[pbc_index].append(gene_cov)
            # add GCs of genes on this contig to all gene GCs
            gene_gc = gene.get_gc_content()
            gene_gcs.append(gene_gc)

    # --------- Compute Contig Stats ------------

    # compute min and max contig length,
    # N50, N90 and total length for all considered contigs
    considered_assembly_length = get_sum_of_contig_lengths(a)
    considered_contigs_n50 = get_Nx(a, 50, contig_lengths)
    considered_contigs_n90 = get_Nx(a, 90, contig_lengths)
    considered_contigs_min_length = min(contig_lengths)
    considered_contigs_max_length = max(contig_lengths)

    # calculate z-score vector for complete assembly (only taking specified contigs into account)
    considered_z_score_vector = \
    calculate_z_score_vector(considered_tetranuc_freqs, considered_trinuc_freqs, considered_dinuc_freqs)

    # compute min and max contig length,
    # N50, N90 and total length for the complete assembly
    total_assembly_length = sum(a.get_all_contig_lengths())
    all_contigs_n50 = get_Nx(a, 50, a.get_all_contig_lengths())
    all_contigs_n90 = get_Nx(a, 90, a.get_all_contig_lengths())
    all_contigs_min_length = min(a.get_all_contig_lengths())
    all_contigs_max_length = max(a.get_all_contig_lengths())

    # compute min and max, weighted mean and SD over all contig covs
    contig_min_cov = {pbc_index: min(coverages)
                        for pbc_index, coverages in contig_coverages.items()}
    contig_max_cov = {pbc_index: max(coverages)
                        for pbc_index, coverages in contig_coverages.items()}
    # WEIGHTING
    all_contigs_mean_cov = {}
    all_contigs_cov_sd = {}
    for pbc_index in a.get_pbc_paths().keys():
        mean_cov, cov_sd = lgth_corr_function( \
                            contig_coverages[pbc_index], \
                            contig_covered_bases_lengths[pbc_index])
        all_contigs_mean_cov[pbc_index] = mean_cov
        all_contigs_cov_sd[pbc_index] = cov_sd

    # compute min and max, weighted mean and SD over all contig GCs
    contig_min_gc = min(contig_gcs)
    contig_max_gc = max(contig_gcs)
    # WEIGHTING
    all_contigs_mean_gc, all_contigs_gc_sd \
    = lgth_corr_function(contig_gcs, contig_lengths)

    # --------- Compute Gene Stats ------------
    # compute mean, SD, min and max over all gene lengths
    all_genes_mean_length = np.mean(gene_lengths)
    all_genes_length_sd = np.std(gene_lengths)
    gene_min_length = min(gene_lengths)
    gene_max_length = max(gene_lengths)

    # compute min and max, weighted mean and SD over all gene covs

    gene_min_cov = {pbc_index: min(coverages)
                        for pbc_index, coverages in gene_coverages.items()}
    gene_max_cov = {pbc_index: max(coverages)
                        for pbc_index, coverages in gene_coverages.items()}
    # WEIGHTING
    all_genes_mean_cov = {}
    all_genes_cov_sd = {}
    for pbc_index in a.get_pbc_paths().keys():
        mean_cov, cov_sd = lgth_corr_function( \
                            gene_coverages[pbc_index], \
                            gene_covered_bases_lengths[pbc_index])
        all_genes_mean_cov[pbc_index] = mean_cov
        all_genes_cov_sd[pbc_index] = cov_sd
    # UNWEIGHTED COV MEAN AND SD
    all_genes_mean_cov_unweighted = {}
    all_genes_cov_sd_unweighted = {}
    for pbc_index in a.get_pbc_paths().keys():
        mean_cov = np.mean(gene_coverages[pbc_index])
        cov_sd = np.std(gene_coverages[pbc_index])
        all_genes_mean_cov_unweighted[pbc_index] = mean_cov
        all_genes_cov_sd_unweighted[pbc_index] = cov_sd

    # compute min and max, weighted man and SD over all gene GCs
    gene_min_gc = min(gene_gcs)
    gene_max_gc = max(gene_gcs)
    # WEIGHTING
    all_genes_mean_gc, all_genes_gc_sd \
    = lgth_corr_function(gene_gcs, gene_lengths)

    # ---------  OUTPUT  ------------

    output_dict = {"contig gc mean" : all_contigs_mean_gc,
        "contig gc sd" : all_contigs_gc_sd,
        "contig gc min" : contig_min_gc,
        "contig gc max" : contig_max_gc,

         "z-score vector (considered contigs)" : considered_z_score_vector,
         "z-score vector (all contigs)" : a.get_total_z_score_vector(),


         "total assembly length" : total_assembly_length,
         "all contigs n50" : all_contigs_n50,
         "all contigs n90" : all_contigs_n90,
         "all contigs len min" : all_contigs_min_length,
         "all contigs len max" : all_contigs_max_length,

        "considered assembly length" : considered_assembly_length,
         "considered contigs n50" : considered_contigs_n50,
         "considered contigs n90" : considered_contigs_n90,
         "considered contigs len min" : considered_contigs_min_length,
         "considered contigs len max" : considered_contigs_max_length,

         "gene gc mean" : all_genes_mean_gc,
         "gene gc sd" : all_genes_gc_sd,
         "gene gc min" : gene_min_gc,
         "gene gc max" : gene_max_gc,

          "gene len mean" : all_genes_mean_length,
          "gene len sd" : all_genes_length_sd,
          "gene len min" : gene_min_length,
          "gene len max" : gene_max_length,

          "contig cov mean" : all_contigs_mean_cov,
          "contig cov sd" : all_contigs_cov_sd,
          "contig cov min" : contig_min_cov,
          "contig cov max" : contig_max_cov,

          "gene cov mean" : all_genes_mean_cov,
          "gene cov sd" : all_genes_cov_sd,
          "gene cov min" : gene_min_cov,
          "gene cov max" : gene_max_cov,
          "gene cov mean unweighted" : all_genes_mean_cov,
          "gene cov sd unweighted" : all_genes_cov_sd,
     }
    # for i in range(len(pbc_paths)):
    #     cov_dict = {"contig cov mean "+str(i) : all_contigs_mean_cov[i],
    #         "contig cov sd "+str(i) : all_contigs_cov_sd[i],
    #          "contig cov min "+str(i) : contig_min_cov[i],
    #          "contig cov max "+str(i) : contig_max_cov[i],
    #
    #          "gene cov mean "+str(i) : all_genes_mean_cov[i],
    #          "gene cov sd "+str(i) : all_genes_cov_sd[i],
    #          "gene cov min "+str(i) : gene_min_cov[i],
    #          "gene cov max "+str(i) : gene_max_cov[i]}
    #     output_dict.update(cov_dict)

    return output_dict

def filter_contigs_and_genes(a, cfg):
    """
    Stores names of all contigs without genes in a list and returns it.
    Moves all contigs without coverage info to a new dict and returns it.
    Moves all genes without coverage info to a new dict and returns it.

    Args:
      a:

    Returns:

    """

    for contig_name, contig in a.get_contigs().items():
        # store contig as geneless if true
        if contig.geneless_info() == 1:
            a.add_geneless_contig(contig_name, contig)


        # if no per base coverage info is available for this contig
        if cfg.include_coverage and contig.no_coverage_info():
            # store it to appropriate dict
            a.add_contig_without_cov(contig_name, contig)

            # and store all its genes to an appropriate dict
            for gene_name in contig.get_genes():
                a.add_gene_without_cov(gene_name, a.get_gene(gene_name))

    if cfg.include_coverage:
        # delete all contigs & genes that have no cov info
        for gene_name, gene in a.get_genes_without_cov().items():
            a.remove_gene(gene_name)

        for contig_name, contig in a.get_contigs_without_cov().items():
            a.remove_contig(contig_name)


def get_raw_array(a, stats_ref):
    """
    Loops over all  contigs requested in a.contigs.keys()
    and returns their metrics in an (n x p) "matrix" (data type: array of arrays)
    where each of the n arrays (rows) describes a gene
    on p variables (columns; i.e. array has p elements).
    output_array contains "raw" data, i.e. prior to NaN imputatioin and rescaling.

    Args:
      a:
      stats_ref:

    Returns:

    """

    output_array = []

    for contig_name, contig in a.get_contigs().items():

        c_ngenes = contig.get_number_of_genes()

        # only if the contig has any genes:
        if c_ngenes > 0:

            c_len = contig.get_length()
            c_pct_assembly_len = percentage_assembly_length(c_len, stats_ref["considered assembly length"])
            # to get the percentage of only the contigs considered, this line should be:
            #c_pct_assembly_len = percentage_assembly_length(c_len, stats_ref["considered assembly length"])

            # -------oligo info--------
            # compute correlation between the z-score vectors of the oligonucleotide frequencies
            # of this contig and of all contigs (only those considered)
            c_z_score_vector = contig.get_z_score_vector()
            c_pearson_test = pearsonr_with_nans_omitted(c_z_score_vector,
                                                  stats_ref["z-score vector (considered contigs)"])
            c_pearson_r = c_pearson_test[0] # get the correlation coefficient
            c_pearson_p = c_pearson_test[1] # get the p-value of the correlation
            # -------oligo info--------

            c_gc_content = contig.get_gc_content()
            c_gcdev = contig.gcdev_from_overall(stats_ref["contig gc mean"],
                                                   stats_ref["contig gc sd"])

            c_cov = {pbc_index: contig.get_coverage(pbc_index) for pbc_index in a.get_pbc_paths().keys()}
            c_covsd = {pbc_index: contig.get_coverage_sd(pbc_index) for pbc_index in a.get_pbc_paths().keys()}
            c_covdev = {pbc_index: contig.covdev_from_overall(stats_ref["contig cov mean"][pbc_index],
                            stats_ref["contig cov sd"][pbc_index], pbc_index)
                            for pbc_index in a.get_pbc_paths().keys()}
            c_genecovm = {pbc_index: contig.get_gene_coverage_mean(pbc_index) for pbc_index in a.get_pbc_paths().keys()}
            c_genecovsd = {pbc_index: contig.get_gene_coverage_sd(pbc_index) for pbc_index in a.get_pbc_paths().keys()}
            c_genelenm = contig.get_gene_lengths_mean()
            c_genelensd = contig.get_gene_lengths_sd()

            for gene_name in contig.get_genes():
                gene = a.get_gene(gene_name)

                # compute correlation between the z-score vectors of the oligonucleotide frequencies ...
                #        ... of this gene and of all contigs (only those considered)
                g_z_score_vector = gene.get_z_score_vector()
                g_pearson_test_o = pearsonr_with_nans_omitted(g_z_score_vector,
                                        stats_ref["z-score vector (considered contigs)"])
                g_pearson_r_o = g_pearson_test_o[0] # get the correlation coefficient
                g_pearson_p_o = g_pearson_test_o[1] # get the p-value of the correlation

                #        ... of this gene and of its contig
                g_pearson_test_c = pearsonr_with_nans_omitted(g_z_score_vector, c_z_score_vector)
                g_pearson_r_c = g_pearson_test_c[0] # get the correlation coefficient
                g_pearson_p_c = g_pearson_test_c[1] # get the p-value of the correlation
                # ------------------

                gene_array = [

                    gene_name, # index 0
                    contig_name, # index 1

                    # spatial contig variables:
                    c_ngenes, # index 2
                    c_len, # index 3
                    c_pct_assembly_len, # c_lendev # index 4
                    c_genelenm, # index 5
                    c_genelensd, # index 6

                    # coverage-related contig variables:
                    c_cov,  # index 7
                    c_covsd, # index 8
                    c_covdev, # index 9
                    c_genecovm, # index 10
                    c_genecovsd, # index 11

                    # compositional contig variables:
                    c_pearson_r, # index 12
                    c_pearson_p, # index 13
                    c_gc_content, # index 14
                    c_gcdev, # index 15

                    # spatial gene variables:
                    gene.get_length(), # index 16
                    gene.lendev_from_contig(a), # index 17
                    gene.lendev_from_overall(stats_ref["gene len mean"],
                                            stats_ref["gene len sd"]), # index 18
                    gene.get_absolute_pos(), # index 19
                    gene.get_terminal_info(), # index 20
                    gene.get_single_gene_info(), # index 21

                    # coverage-related gene variables:
                    {pbc_index: gene.get_coverage(pbc_index) for pbc_index in a.get_pbc_paths().keys()}, # index 22
                    {pbc_index: gene.get_coverage_sd(pbc_index) for pbc_index in a.get_pbc_paths().keys()}, # index 23
                    {pbc_index: gene.covdev_from_contig(a, pbc_index) for pbc_index in a.get_pbc_paths().keys()}, # index 24
                    {pbc_index: gene.covdev_from_overall(stats_ref["gene cov mean"][pbc_index], stats_ref["gene cov sd"][pbc_index], pbc_index) \
                        for pbc_index in a.get_pbc_paths().keys()}, # index 25

                    # compositional gene variables:
                    g_pearson_r_o, # index 26
                    g_pearson_p_o, # index 27
                    g_pearson_r_c, # index 28
                    g_pearson_p_c, # index 29

                    gene.get_gc_content(), # index 30
                    gene.gcdev_from_contig(a), # index 31
                    gene.gcdev_from_overall(stats_ref["gene gc mean"],
                                            stats_ref["gene gc sd"]), # index 32
                ]

                output_array.append(gene_array)

    return output_array


def output_table(a, raw_array, filename):
    """
    Prints raw_array including header to a file.

    Args:
      a:
      raw_array:
      filename:

    Returns:

    """

    path1_raw_array = a.get_output_dir() + "/" + filename + ".txt"
    path2_raw_array = a.get_output_dir() + "/" + filename + ".csv"
    out1_raw_array = open(path1_raw_array, "w")
    out2_raw_array = open(path2_raw_array, "w")

    c_cov_header = []
    g_cov_header = []
    for elem in ["c_cov", "c_covsd", "c_covdev", "c_genecovm", "c_genecovsd"]:
        for pbc_index in sorted(a.get_pbc_paths().keys()):
            c_cov_header.append(elem+"_"+str(pbc_index))
    for elem in (["g_cov", "g_covsd", "g_covdev_c", "g_covdev_o"]):
        for pbc_index in sorted(a.get_pbc_paths().keys()):
                g_cov_header.append(elem+"_"+str(pbc_index))

    header = [
        "g_name", "c_name",
        "c_num_of_genes", "c_len", "c_pct_assemby_len", "c_genelenm", "c_genelensd"] + \
        c_cov_header + \
        ["c_pearson_r", "c_pearson_p", "c_gc_cont", "c_gcdev",
        "g_len", "g_lendev_c", "g_lendev_o", "g_abspos", "g_terminal", "g_single"] + \
        g_cov_header + \
        ["g_pearson_r_o", "g_pearson_p_o", "g_pearson_r_c", "g_pearson_p_c",
        "g_gc_cont", "g_gcdev_c", "g_gcdev_o" # , "g_single_exon"
    ]

    out1_raw_array.write(("\t".join(header)) + "\n")
    out2_raw_array.write((",".join(header)) + "\n")

    for gene_array in raw_array:
        flat_gene_array = []
        for item in gene_array:
            if type(item) == dict:
                for subitem in sorted(item.keys()):
                    flat_gene_array.append(item.get(subitem))
            else:
                flat_gene_array.append(item)
        if len(header) != len(flat_gene_array):
            print(header)
            print(len(header), len(flat_gene_array))
            print(flat_gene_array)
            print(gene_array)
        out1_raw_array.write(("\t".join(str(round_if_float(x, 4)) for x in flat_gene_array)) + "\n")
        out2_raw_array.write((",".join(str(round_if_float(x, 4)) for x in flat_gene_array)) + "\n")

    out1_raw_array.close()
    out2_raw_array.close()


def complementary_assembly_length(c_len, total_assembly_length):
    """
    Returns how much percent of the total assembly length is
    NOT made up by this contig.
    E.g.   total length 100, contig length 10 --> will return 0.9
           total length 100, contig length 95 --> will return 0.05
    By this, larger contigs will be penalised less than smaller ones.

    Args:
      c_len:
      total_assembly_length:

    Returns:

    """

    pct_of_total = c_len / total_assembly_length
    inv_pct_of_total = 1 - pct_of_total
    return inv_pct_of_total


def percentage_assembly_length(c_len, total_assembly_length):
    """
    Returns the percentage of total assembly length
    covered by this contig.

    Args:
      c_len:
      total_assembly_length:

    Returns:

    """

    pct_of_total = c_len / total_assembly_length
    return pct_of_total


def deviation_from_n50(c_len, n50):
    """
    Returns the contig's deviation from N50,
    normalised by N50.
    If the contig is smaller than N50 this will be a positive value between 0 and 1.
    If the contig is larger than N50 this will be a negative value.
    By adding the value (NOT the absolute value) to the sum of metrics,
    smaller contigs are penalised while larger contigs are rewarded.

    Args:
      c_len:
      n50:

    Returns:

    """

    deviation = n50 - c_len
    normalised_deviation = deviation / n50
    return normalised_deviation


def output_summary_file(a, stats_ref, filename, include_coverage):
    """

    Args:
      a:
      stats_ref:
      filename:
      include_coverage:

    Returns:

    """
    path_summaryfile = a.get_output_dir() + "/" + filename
    summaryfile = open(path_summaryfile, "w")

    contig_cov_summary = ""
    gene_cov_summary = ""
    if include_coverage:
        for pbc_index in a.get_pbc_paths().keys():
            contig_cov_summary += "contig cov mean\t" + str(round(stats_ref["contig cov mean"][pbc_index], 2)) + "\n" + \
             "contig cov sd\t" + str(round(stats_ref["contig cov sd"][pbc_index], 2)) + "\n" + \
             "contig cov min\t" + str(round(stats_ref["contig cov min"][pbc_index], 2)) + "\n" + \
             "contig cov max\t" + str(round(stats_ref["contig cov max"][pbc_index], 2)) + "\n"
        contig_cov_summary += "\n"

        for pbc_index in a.get_pbc_paths().keys():
            gene_cov_summary += "gene cov mean\t" + str(round(stats_ref["gene cov mean"][pbc_index], 2))  + "\n" + \
              "gene cov sd\t" + str(round(stats_ref["gene cov sd"][pbc_index], 2)) + "\n" + \
              "gene cov min\t" + str(round(stats_ref["gene cov min"][pbc_index], 2))  + "\n" + \
              "gene cov max\t" + str(round(stats_ref["gene cov max"][pbc_index], 2)) + "\n"
        gene_cov_summary += "\n"

    summaryfile.write("total # of contigs:\t" + str(len(a.get_contigs()) + len(a.get_contigs_without_cov())) + "\n" +
                      "\t# geneless contigs:\t" + str(len(a.get_geneless_contigs()))  + "\n" +
                      "\t# contigs w/ genes:\t"
                          + str(len(a.get_contigs()) + len(a.get_contigs_without_cov()) - len(a.get_geneless_contigs())) + "\n" +

                      "\t# contigs w/o cov info:\t" + str(len(a.get_contigs_without_cov())) + "\n" +
                      "\t# contigs w/ cov info\t" + str(len(a.get_contigs())) + "\n" +

                    "\t# geneless contigs w/o cov info:\t"
                                  + str(len(set(a.get_contigs_without_cov()).intersection(a.get_geneless_contigs()))) + "\n" +
                    "\t# contigs with genes and w/ cov info:\t"
                              + str(len(a.get_contigs()) - len(a.get_geneless_contigs())) + "\n\n" +

                "total assembly length\t" + str(stats_ref["total assembly length"]) + "\n" +
                 "all contigs N50\t" + str(round(stats_ref["all contigs n50"], 2)) + "\n" +
                  "all contigs N90\t" + str(round(stats_ref["all contigs n90"], 2)) + "\n" +
                  "all contigs len min\t" + str(stats_ref["all contigs len min"]) + "\n" +
                  "all contigs len max\t" + str(stats_ref["all contigs len max"]) + "\n\n" +

                "considered assembly length\t" + str(stats_ref["considered assembly length"]) + "\n" +
                "considered contigs N50\t" + str(round(stats_ref["considered contigs n50"], 2)) + "\n" +
                "considered contigs N90\t" + str(round(stats_ref["considered contigs n90"], 2)) + "\n" +
                "considered contigs len min\t" + str(stats_ref["considered contigs len min"]) + "\n" +
                "considered contigs len max\t" + str(stats_ref["considered contigs len max"]) + "\n\n" +

                contig_cov_summary +

                "contig gc mean\t" + str(round(stats_ref["contig gc mean"], 2)) + "\n" +
                 "contig gc sd\t" + str(round(stats_ref["contig gc sd"], 2)) + "\n" +
                 "contig gc min\t" + str(round(stats_ref["contig gc min"], 2)) + "\n" +
                 "contig gc max\t" + str(round(stats_ref["contig gc max"], 2)) + "\n\n" +

                "total # genes:\t" + str(len(a.get_genes()) + len(a.get_genes_without_cov())) + "\n" +
                "\t# genes w/o cov info:\t" + str(len(a.get_genes_without_cov())) + "\n" +
                "\t# genes w/ cov info\t" + str(len(a.get_genes()))  + "\n\n" +

                 "gene len mean\t" + str(round(stats_ref["gene len mean"], 2))  + "\n" +
                  "gene len sd\t" + str(round(stats_ref["gene len sd"], 2)) + "\n" +
                  "gene len min\t" + str(stats_ref["gene len min"]) + "\n" +
                  "gene len max\t" + str(stats_ref["gene len max"]) + "\n\n" +

                  gene_cov_summary +

                  "gene gc mean\t" + str(round(stats_ref["gene gc mean"], 2))  + "\n" +
                  "gene gc sd\t" + str(round(stats_ref["gene gc sd"], 2)) + "\n" +
                  "gene gc min\t" + str(round(stats_ref["gene gc min"], 2))  + "\n" +
                  "gene gc max\t" + str(round(stats_ref["gene gc max"], 2)) + "\n"
                )

    summaryfile.close()


# !!!!!!!!!!!!!!!!!!!!!!!!!!! READING INPUT FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!


# ==============================================================
# ======================= READ GFF FILE ========================
# ==============================================================

def init_contig(a, name, length):
    """
    initalize contig with name and length

    Args:
      a:
      name:
      length:
    """
    # add this length to a.all_contig_lengths
    # to store the size of the initial (unfiltered) assembly
    a.add_contig_length(length)
    # initialise contig
    contig = Contig(name, length, a)
    # add to dictionary of contigs
    a.add_contig(name, contig)


def read_faidx(a):
    """ read FASTA index to initalize contigs with length information provided

    Args:
      a:
    """

    with open(a.get_output_path()+"tmp/tmp.MILTS.fasta.fai", 'r') as faidx: # open FAIDX file
        for line in faidx:
            line_array = line.split()
            contig_name = line_array[0]
            if contig_name in a.get_contigs().keys():
                continue
            else:
                contig_length = int(line_array[1])
                init_contig(a, contig_name, contig_length)


def read_gff(a, gene_tag, source_type, include_pseudogenes):
    """

    Args:
      a:
      gene_tag:
      source_type:
      include_pseudogenes:

    Returns:

    """

    with open(a.get_gff_path(), 'r') as gff: # open GFF file
        for line in gff:
            if not line.startswith('#'): # if the line holds annotations
                tmp_array = line.split()

                if ((source_type == "default" and tmp_array[2] == gene_tag) or \
                   (tmp_array[1] == source_type and tmp_array[2] == gene_tag) or \
                   (tmp_array[2] == "pseudogene" and include_pseudogenes == True)):



                    associated_contig = tmp_array[0]
                    source = tmp_array[1]
                    start_pos = int(tmp_array[3])
                    end_pos = int(tmp_array[4])

                    if tmp_array[5] == '.':          # if there is no score
                        score = tmp_array[5]         # simply add the '.'
                    else:                            # if a score is given
                        score = float(tmp_array[5])  # add it as a float

                    strand = tmp_array[6]
                    attributes = tmp_array[8].split(';')
                    id_attribute = attributes[0].split('=') # something like ID=gene1
                    gene_name = strip_ID(id_attribute[1]) # get only the name after the equals sign without "gene:" or "gene-" prefix

                    # initialise gene
                    gene = Gene(gene_name, start_pos, end_pos, associated_contig, source, score, strand, attributes, a)

                    if 'part=' in tmp_array[8]:
                        # dont add partial genes
                        # TODO: add ability to handle partial genes
                        a.add_partial_gene(gene_name, gene)
                        continue

                    # add to dictionary of genes
                    a.add_gene(gene_name, gene)

                    # add to list of genes in associated contig
                    a.get_contig(associated_contig).add_gene(gene_name)
                else:
                    pass
            elif "#FASTA" in line: # if FASTA block has been reached
                break            # stop parsing


# ==============================================================
# ===================== READ FASTA FILE =====================
# ==============================================================


# ---- start: variables needed for oligonucleotide frequencies ------
def init_tetranuc_freq():
    """ """

    # intitialise dictionaries to keep count of overall oligo frequencies
    total_tetranuc_freqs = {nuc : 0 for nuc in TETRANUCLEOTIDES}
    total_trinuc_freqs = {nuc : 0 for nuc in TRINUCLEOTIDES}
    total_dinuc_freqs = {nuc : 0 for nuc in DINUCLEOTIDES}

    return total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs
# ---- end: variables needed for oligonucleotide frequencies ------


def return_positions_of_Ns(sequence):
    """For a given sequence (e.g. scaffold / contig / gene) this function
    will return a set holding all indices (1-)

    Args:
      sequence:

    Returns:

    """
    return {(i+1) for i, base in enumerate(sequence) if base == "N"}


def set_seq_info(a, current_contig, raw_seq, total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs):
    """set sequence info for given contig.

    Args:
      a:
      current_contig:
      raw_seq:
      total_tetranuc_freqs:
      total_trinuc_freqs:
      total_dinuc_freqs:

    Returns:

    """
    # get a set holding the (1-based) positions of all Ns in this contig
    contig_N_pos = return_positions_of_Ns(raw_seq)
    # pass the set to the contig
    current_contig.set_positions_of_Ns(contig_N_pos)
    contig_seq = raw_seq.replace("N", "") # remove all ambiguous characters
    # CAUTION!! LENGTH OF CONTIGS / GENES IS *WITH* AMBIGUOUS CHARACTERS

    # -----> OLIGOTEST
    # count all tetra-/tri-/dinucleotides in the contig sequence and assign to contig
    contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs = \
      get_all_obsfreqs(contig_seq)
    current_contig.set_oligofreqs(contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)
    # calculate and set the vector of the 256 oligonucleotide freq z-scores of this contig
    contig_z_score_vector =\
      calculate_z_score_vector(contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)
    current_contig.set_z_score_vector(contig_z_score_vector)
    # add oligofrequencies in contig to total oligofrequencies observed so far
    total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs = \
      update_observed_oligofreqs(total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs,
                                    contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)
    # <----- OLIGOTEST

    # compute GC content of this contig
    contig_gc = (contig_seq.count('G')  / len(contig_seq)) \
      + (contig_seq.count('C')  / len(contig_seq))

    current_contig.set_gc_content(contig_gc)

    # for each gene on this contig
    for gene_name in current_contig.get_genes():
        current_gene = a.get_gene(gene_name)

        # ---- get start and end position ------------
        # because bp are 1-based but python lists are 0-based:
        gene_seq_start = current_gene.get_start_pos() - 1
        # because list slicing ([cov_start:cov_end]) stops BEFORE the given index:
        gene_seq_end = current_gene.get_end_pos()
        # (one would want to have gene_end+1, but since the array is 0-based
        # it's already "shifted" +1 to the right
        # --------------------------------------------

        # regard only gene region on contig
        gene_seq = raw_seq[gene_seq_start:gene_seq_end].replace("N", "")

        # count all tetra-/tri-/dinucleotides in the gene sequence and assign to gene
        gene_tetranuc_freqs, gene_trinuc_freqs, gene_dinuc_freqs = \
          get_all_obsfreqs(gene_seq)
        current_gene.set_oligofreqs(gene_tetranuc_freqs, gene_trinuc_freqs, gene_dinuc_freqs)

        # calculate and set the vector of the 256 oligonucleotide freq z-scores of this gene
        gene_z_score_vector =\
          calculate_z_score_vector(gene_tetranuc_freqs, gene_trinuc_freqs, gene_dinuc_freqs)
        current_gene.set_z_score_vector(gene_z_score_vector)

        # compute GC content of current gene
        gene_gc = (gene_seq.count('G') / len(gene_seq)) \
          + (gene_seq.count('C') / len(gene_seq))

        current_gene.set_gc_content(gene_gc)

    return total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs


def read_fasta(a):
    """

    Args:
      a:

    Returns:

    """

    # variable needed to keep track of currently parsed contig
    current_contig = None
    total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs = init_tetranuc_freq()

    with open(a.get_fasta_path(), 'r') as fasta: # open FASTA file
        for line in fasta:
            line = line.strip()
            # if line contains fasta header
            # store sequence of previous contig and save new name
            if line.startswith(">"):

                if current_contig:
                    total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs = \
                      set_seq_info(a, current_contig, raw_seq, total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs)

                # refresh for new contig
                raw_seq = ''
                tmp_array = line.split()
                tmp_contig_name = tmp_array[0][1:]
                current_contig = a.get_contig(tmp_contig_name)
                if current_contig.geneless_info() == 1: # if contig is geneless
                    # a.remove_contig(tmp_contig_name)     #TODO: remove from list?!?!
                    a.add_geneless_contig(tmp_contig_name, current_contig)
                    current_contig = None               # don't process further

            elif current_contig: # current contig is assigned
                # just in case (for case-insensitive counting)
                raw_seq = raw_seq + line.upper() # convert all bases to upper-case

            else:
                # do nothing if current contig is not annotated in GFF (geneless)
                continue

    if current_contig:
        # handle last contig
        total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs = \
          set_seq_info(a, current_contig, raw_seq, total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs)

    a.set_total_z_score_vector(calculate_z_score_vector(total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs))


# ==============================================================
# ===================== READ COVERAGE FILE =====================
# ==============================================================

def read_pbc(a, include_coverage):
    """

    Args:
      a:
      include_coverage:

    Returns:

    """

    current_contig = None # contig (object) to which coverages are currently added
    current_contig_name = "dummy" # name of contig to which coverages are currently added

    if include_coverage:
        for pbc_index, pbc_path in a.get_pbc_paths().items():
            with open(pbc_path, 'r') as pbc: # open PBC file
                for line in pbc:
                    pbc_info_array = line.split()

                    if not pbc_info_array[0] == current_contig_name: # if new contig has been reached
                        current_contig_name = pbc_info_array[0]      # update name
                        if current_contig_name in a.get_contigs().keys():       # if this is a contig with genes
                            current_contig = a.get_contig(current_contig_name) # update current contig
                        else: # if current contig is not in GFF contig list (i.e. has no genes)
                            current_contig = None # set this to None so that all other info for this contig gets ignored
                            """print("No coverage information for contig ", current_contig_name,
                                  " stored since it has no genes.", file=open("error.txt", "a"))"""
                    if current_contig: # only add to contigs which exist in GFF file
                        base = int(pbc_info_array[1])
                        coverage = int(float(pbc_info_array[2]))
                        current_contig.add_base_coverage(pbc_index, base, coverage)
    else: # if user does not want to include coverage information or it is not available
        mock_coverage = 1 # coverage that is set for all bases
        for contig in a.get_contigs().values():
            if contig.geneless_info() == 0: # only store "info" if contig has genes
                for base in range(1, contig.get_length()+1): # set coverage for every base to mock coverage
                    for pbc_index in a.get_pbc_paths().keys():
                        contig.add_base_coverage(pbc_index, base, mock_coverage)



# # !!!!!!!!!!!!!!!!!!!!!!!!!!! GENERATE OUTPUT PART !!!!!!!!!!!!!!!!!!!!!!!!!!!

def process_gene_info(cfg):
    """

    Args:
      cfg:

    Returns:

    """

    # TODO: add to logging info in beginning
    # logging.debug("include pseudogenes = " + str(cfg.include_pseudogenes))
    # logging.debug("include coverage = " + str(cfg.include_coverage))

    # ====================== VARIABLES ======================

    a = Assembly(cfg.gff_path, cfg.fasta_path, cfg.pbc_paths, cfg.output_path)

    read_faidx(a) # use faidx to initalize contigs
    read_gff(a, cfg.gff_gene_tag, cfg.gff_source, cfg.include_pseudogenes)
    if a.partial_genes:
        logging.warning('Gene(s) excluded due to partialness:\n{}'.format(list(a.partial_genes.keys())))
    read_fasta(a)
    read_pbc(a, cfg.include_coverage)

    # store names of geneless contigs to a list
    # exclude all contigs and genes without coverage from further processing
    # (i.e. remove them from a.contigs and store them to respective own dicts)
    filter_contigs_and_genes(a, cfg)
    # AT THIS POINT, a.contigs WILL ONLY CONTAIN CONTIGS WITH COVERAGE INFO
    # (BUT MAY STILL CONTAIN GENELESS CONTIGS)

    # use per base coverage information to compute coverage mean & SD for all genes and contigs
    # also compute absolute positions of genes on contigs
    compute_coverage_and_positional_info(a)

    # ----------------------- GET STATS -> values for reference ----------------------

    # get stats taking ALL contigs and genes into account (except those without coverage)
    stats_for_all_contigs = compute_stats(a, get_lgth_corr_stats)

    # # get stats taking all contigs & genes into account but WEIGHTING BY SQUARED LENGTH
    # stats_sq_weighting = compute_stats(a.contigs.keys(), get_lgth_corr_stats_weighted_by_sq_lgth)

    # # long contigs (> 1000 or >600, resp.) only --------------------------
    # contigs1000 = []
    # contig1000_genes = []
    # for contig_name, contig in a.contigs.items():
    #     if contig.get_length() > 1000:
    #         contigs1000.append(contig_name)
    #         contig1000_genes.extend(contig.get_genes())

    # stats_for_contigs1000 = compute_stats(contigs1000, get_lgth_corr_stats)

    # ----------------------------------------------------------------------------------

    pathlib.Path(a.get_output_dir()).mkdir(parents=True, exist_ok=True)

    raw_array = get_raw_array(a, stats_for_all_contigs)
    output_table(a, raw_array, "raw_gene_table")
    imputed_array = impute_array(a, raw_array, cfg.include_coverage)
    output_table(a, imputed_array, "imputed_gene_table")

    output_summary_file(a, stats_for_all_contigs, "summary.txt", cfg.include_coverage)

def main():
    """ """
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(config_path)

    process_gene_info(cfg)

if __name__ == '__main__':
    main()
