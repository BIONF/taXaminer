#!/usr/bin/env python

"""Prepare config parameters and check data validity.

Processes the users config file, sets default values for those not
given by the user, checks file existences to assess what needs to be
computed and makes sanity checks for input data

Expects path to config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

import pathlib # to create directories
import yaml # read config file
import sys # parse command line arguments
import logging
import taxopy



class Config:
    """Object to hold all configuration options
    """

    def __init__(self, cfg_dict):
        """
        Constructs configuration object.

        Parameters
        ----------
            cfg_dict : dict
                dictionary holding the preprocessed configuration parameters
        """
        self.fasta_path = cfg_dict.get('fasta_path')
        self.gff_path = cfg_dict.get('gff_path')
        self.output_path = cfg_dict.get('output_path')
        self.taxon_id = cfg_dict.get('taxon_id')
        self.taxon_id_rank = cfg_dict.get('taxon_id_rank')
        self.taxon_name = cfg_dict.get('taxon_name')

        self.gff_ta_path = cfg_dict.get('gff_ta_path')
        self.threads = cfg_dict.get('threads')

        self.read_paths = cfg_dict.get('read_paths')
        self.bam_paths = cfg_dict.get('bam_paths')
        self.pbc_paths = cfg_dict.get('pbc_paths')
        self.cov_set_exists = cfg_dict.get('cov_set_exists')
        self.include_coverage = cfg_dict.get('include_coverage')
        self.compute_coverage = cfg_dict.get('compute_coverage')
        self.min_insert = cfg_dict.get('min_insert')
        self.max_insert = cfg_dict.get('max_insert')
        self.read_orientation = cfg_dict.get('read_orientation')

        self.proteins_path = cfg_dict.get('proteins_path')
        self.extract_proteins = cfg_dict.get('extract_proteins')
        self.use_phase_info = cfg_dict.get('use_phase_info')
        self.assignment_mode = cfg_dict.get('assignment_mode')
        self.quick_mode_search_rank = cfg_dict.get('quick_mode_search_rank')
        self.quick_mode_match_rank = cfg_dict.get('quick_mode_match_rank')
        self.tax_assignment_path = cfg_dict.get('tax_assignment_path')
        self.compute_tax_assignment = cfg_dict.get('compute_tax_assignment')
        self.db_dir = f"{cfg_dict.get('db_dir')}/"
        self.database_path = cfg_dict.get('database_path')
        self.taxon_exclude = cfg_dict.get('taxon_exclude')
        self.exclusion_rank = cfg_dict.get('exclusion_rank')

        self.update_plots = cfg_dict.get('update_plots')
        self.num_groups_plot = cfg_dict.get('num_groups_plot')
        self.merging_labels = cfg_dict.get('merging_labels')
        self.output_pdf = cfg_dict.get('output_pdf')
        self.output_png = cfg_dict.get('output_png')

        self.include_pseudogenes = cfg_dict.get('include_pseudogenes')
        self.gff_preset = cfg_dict.get('gff_preset')
        self.input_variables = cfg_dict.get('input_variables')
        self.perform_parallel_analysis = cfg_dict.get('perform_parallel_analysis')
        self.num_pcs = cfg_dict.get('num_pcs')
        self.coverage_cutoff_mode = cfg_dict.get('coverage_cutoff_mode')

        self.perform_kmeans = cfg_dict.get('perform_kmeans')
        self.kmeans_k = cfg_dict.get('kmeans_k')
        self.perform_hclust = cfg_dict.get('perform_hclust')
        self.hclust_k = cfg_dict.get('hclust_k')
        self.perform_mclust = cfg_dict.get('perform_mclust')
        self.mclust_k = cfg_dict.get('mclust_k')
        self.perform_dbscan = cfg_dict.get('perform_dbscan')
        self.dbscan_groups = cfg_dict.get('dbscan_groups')
        self.custom_eps = cfg_dict.get('custom_eps')
        self.custom_minPts = cfg_dict.get('custom_minPts')

        self.script_dir = cfg_dict.get('script_dir')
        self.usr_cfg_path = cfg_dict.get('usr_cfg_path')
        self.cfg_path = cfg_dict.get('cfg_path')

        # GFF parsing rule
        self.gff_source = cfg_dict.get('gff_source')
        self.gff_gene_type = cfg_dict.get('gff_gene_type')
        self.gff_transcript_type = cfg_dict.get('gff_transcript_type') \
                            if 'gff_transcript_type' in cfg_dict.keys() \
                            else cfg_dict.get('gff_fasta_header_type')
        self.gff_cds_type = cfg_dict.get('gff_cds_type')
        self.gff_fasta_header_type = cfg_dict.get('gff_fasta_header_type') \
                            if 'gff_fasta_header_type' in cfg_dict.keys() \
                            else cfg_dict.get('gff_transcript_type')
        self.gff_fasta_header_attr = cfg_dict.get('gff_fasta_header_attr')
        self.gff_connection = cfg_dict.get('gff_connection')
        self.gff_parent_child_types = cfg_dict.get('gff_parent_child_types')
        self.gff_parent_attr = cfg_dict.get('gff_parent_attr')
        self.gff_child_attr = cfg_dict.get('gff_child_attr')
        self.gff_gene_attr = cfg_dict.get('gff_gene_attr')

        # tool cmds
        self.diamond = cfg_dict.get('diamond')
        self.samtools = cfg_dict.get('samtools')
        self.bedtools = cfg_dict.get('bedtools')
        self.bowtie2 = cfg_dict.get('bowtie2')
        self.krona = cfg_dict.get('krona')


def check_pbc_ids(config_obj):
    """Check if contig IDs in PBC match IDs in GFF.

    Args:
      config_obj(obj): Config class object of config parameters
    """

    rm_pbcs = []

    pbc_ids = set()
    gff_ids = set()

    for pbc_index, pbc_path in config_obj.pbc_paths.items():
        with open(pbc_path, 'r') as pbc:
            # collect IDs of scaffolds from FASTA
            for line in pbc:
                pbc_ids.add(line.split()[0])

        with open(config_obj.gff_path, 'r') as gff:
            # collect IDs of scaffolds from GFF
            for line in gff:
                if not line.startswith('#'):
                    gff_ids.add(line.split()[0])

        # match the two sets of IDs and save which could not be matched
        missing = set()
        for p_id in pbc_ids:
            if not p_id in gff_ids:
                missing.add(p_id)

        # geneless scaffolds may not be annotated in GFF
        # thus check if sets of IDs in FASTA and missing are unequal
        if len(pbc_ids) == len(missing):
            # if the sets are equally long, the IDs could not be matched at all
            logging.warning(
                f'PBC headers and scaffold IDs in GFF do not match for PBC '
                f'file {pbc_index}. Removing from further computations')
            rm_pbcs.append(pbc_index)

    for pbc_index in rm_pbcs:
        del config_obj.pbc_paths[pbc_index]
    if not config_obj.pbc_paths.keys():
        logging.error('No PBC headers and scaffold IDs in GFF match. Please '
                      'recheck your files, restart without stating coverage '
                      'information or set "include_coverage" to "FALSE".')
        sys.exit()


def check_assembly_ids(config_obj):
    """Check if headers in FASTA of assembly match IDs in GFF.

    Args:
      config_obj(obj): Config class object of config parameters
    """

    fasta_ids = set()
    gff_ids = set()

    with open(config_obj.output_path+'tmp/tmp.taxaminer.fasta.fai', 'r') as gfai:
        # collect IDs of scaffolds from FASTA
        for line in gfai:
            fasta_ids.add(line.split()[0])

    with open(config_obj.gff_path, 'r') as gff:
        # collect IDs of scaffolds from GFF
        for line in gff:
            if not line.startswith('#'):
                gff_ids.add(line.split()[0])

    # match the two sets of IDs and save which could not be matched
    missing = set()
    for f_id in fasta_ids:
        if not f_id in gff_ids:
            missing.add(f_id)

    # geneless scaffolds may not be annotated in GFF
    # thus check if sets of IDs in FASTA and missing are unequal
    if len(fasta_ids) == len(missing):
        # if the sets are equally long, the IDs could not be matched at all
        sys.exit('Error: FASTA headers and scaffold IDs in GFF do not match')


def set_gff_parsing_rules(config_obj):
    """Parsing which features in GFF to read.

    Set parsing information based on default presets or user input.
    Defines which features are considered as genes, which as
    transcript, where header for FASTA can be found and how gene and
    transcript can be linked

    Args:
      gff_source(str): either name of a preset or path to file with parsing rule

    Returns:
      (dict) dictionary with parameters for parsing the GFF
    """

    gff_dict = {}

    ## presets for parsing of GFF
    if config_obj.get('gff_preset') == "default" or (config_obj.get('gff_preset') == None and config_obj.get('gff_source') == None):
        gff_dict['gff_source'] = 'default'
        gff_dict['gff_gene_type'] = 'gene'
        gff_dict['gff_transcript_type'] = 'mRNA'
        gff_dict['gff_cds_type'] = 'CDS'
        gff_dict['gff_fasta_header_type'] = 'mRNA,CDS'
        gff_dict['gff_fasta_header_attr'] = 'ID'
        gff_dict['gff_connection'] = 'parent/child'
        gff_dict['gff_parent_child_types'] = 'mRNA,CDS'
        gff_dict['gff_parent_attr'] = 'Parent'
        gff_dict['gff_child_attr'] =  'ID'
    elif config_obj.get('gff_preset') == "maker":
        gff_dict['gff_source'] = 'maker'
        gff_dict['gff_gene_type'] = 'gene'
        gff_dict['gff_transcript_type'] = 'mRNA'
        gff_dict['gff_cds_type'] = 'CDS'
        gff_dict['gff_fasta_header_type'] = 'mRNA'
        gff_dict['gff_fasta_header_attr'] = 'ID'
        gff_dict['gff_connection'] = 'parent/child'
        gff_dict['gff_parent_child_types'] = 'mRNA'
        gff_dict['gff_parent_attr'] = 'Parent'
        gff_dict['gff_child_attr'] =  'ID'
    elif config_obj.get('gff_preset') == "augustus_masked":
        gff_dict['gff_source'] = 'augustus_masked'
        gff_dict['gff_gene_type'] = 'match'
        gff_dict['gff_transcript_type'] = 'mRNA'
        gff_dict['gff_cds_type'] = 'match_part'
        gff_dict['gff_fasta_header_type'] = 'match'
        gff_dict['gff_fasta_header_attr'] = 'Name'
        gff_dict['gff_connection'] = 'inline'
        gff_dict['gff_gene_attr'] = 'ID'
    # read file with information on how to parse GFF
    else:
        gff_dict['gff_source'] = config_obj.get('gff_source')
        gff_dict['gff_gene_type'] = config_obj.get('gff_gene_type')
        gff_dict['gff_transcript_type'] = config_obj.get('gff_transcript_type')
        gff_dict['gff_cds_type'] = config_obj.get('gff_cds_type')
        gff_dict['gff_fasta_header_type'] = config_obj.get('gff_fasta_header_type')
        gff_dict['gff_fasta_header_attr'] = config_obj.get('gff_fasta_header_attr')
        gff_dict['gff_connection'] = config_obj.get('gff_connection')
        gff_dict['gff_parent_child_types'] = config_obj.get('gff_parent_child_types')
        gff_dict['gff_parent_attr'] = config_obj.get('gff_parent_attr')
        gff_dict['gff_child_attr'] = config_obj.get('gff_child_attr')
        gff_dict['gff_gene_attr'] = config_obj.get('gff_gene_attr')

    return gff_dict


def enumerated_key(config_obj, key_name, pre_keys, *default):
    """Get matches for key_name in user config file.

    Retrieve all matches for the enumerated key_name parameters (they
    have a numbered suffix) and set the paths, either as given by user
    or default values. pre_keys is a list of indexes that were
    identified for previous key_name and need to get default values if
    not defined by the user [Is used for retrieving the coverage
    information, thus if there is a set of read data with suffix 1 and
    the user does not define paths for BAM and PBC, default locations
    for these need to be set]

    Args:
      config_obj(dict): dictionary of config parameters from user file
      key_name(str): name of the parameter to be matched
      pre_keys(list): keys for which values need to be defined
      *default(str): default values for the key

    Returns:
      (dict) {enumerated suffix (index) : value}
    """

    # find matches for key
    matches = []
    for key in config_obj.keys():
        if key_name in key:
            matches.append(key)

    # check if there is a digit in the matches
    # to identify the index of the coverage set
    dict = {}
    for match in matches:
        mod_default = None
        if config_obj.get(match):
            # catch user stating e.g. 'pbc_path' as config parameter
            if len(match.split('_')) >= 3 and match.split('_')[2].isdigit():
                match_num = int(match.split('_')[2])
                if type(config_obj[match]) == list:
                    if any(path.startswith('path/to/') for path in config_obj[match]):
                        logging.warning(f'{match}: path can not start with '
                                    f'"path/to/". Trying to use default instead')
                elif type(config_obj[match]) == str:
                    if config_obj[match].startswith('path/to/'):
                        logging.warning(f'{match}: path can not start with '
                                    f'"path/to/". Trying to use default instead')
                if default and '/' in default[0]:
                    mod_default = '.'.join(default[0].split('.')[:-1])+'_'+str(match_num)+'.'+default[0].split('.')[-1]
                dict[match_num] = set_default(config_obj, match,
                                              default if not mod_default else mod_default)
            else:
                dict[0] = set_default(config_obj, match, default)

    for pre_key in pre_keys:
        # if there are more values required than given,
        # fill them up with default paths
        if not dict.get(pre_key):
            # if the default value is a path ('/' in there) then alter the
            # path to append the index of the key_name to it
            if '/' in default[0]:
                dict[pre_key] = '.'.join(default[0].split('.')[:-1])+'_'+str(pre_key)+'.'+default[0].split('.')[-1]
            else:
                dict[pre_key] = default[0]

    return dict


def check_dir_slash(path):
    """Appends trailing slash to path if necessary.

    Args:
      path(str): path to be checked for slash

    Returns:
      (str) path with slash at the end

    """
    if path.endswith('/'):
        return path
    else:
        return path+'/'


def check_file_inexistence(paths, and_or):
    """Check if file(s) at (list of) path(s) do not exists.

    Checks for a list of paths or a single one if these files exist. Returns
    'FALSE' if they do exist, else 'TRUE'. The parameter and_or decides if
    all ('AND') or at least one ('OR') of the need to exist

    Args:
      paths(str,list): path or list of paths to check file existence
      and_or(str): 'AND' if all files must exist, 'OR' if only one has to

    Returns:
        (str) 'FALSE' if the required files do exist, else 'TRUE'
    """

    if type(paths) != list:
        paths = [paths]

    counter_e = 0 # counter of how many files exist
    for path in paths:
        file = pathlib.Path(path)
        if file.is_file():
            counter_e += 1

    # check if all or at least one of the files exist
    if (and_or == 'AND' and counter_e == len(paths)) or \
     (and_or == 'OR' and counter_e > 0):
        # returns 'FALSE' if they do as this means no computations
        # need to be performed
        return 'FALSE'
    else:
        return 'TRUE'


def set_default(config_obj, key_name, *default):
    """Check if parameter is given in config file, else set to default.

    Args:
      config_obj(dict): dict of user input for config parameters
      key_name(str): key to look for in config_obj
      *default(str): default parameter for key_name

    Returns:
        (str) if key_name in config_obj then value, else default

    """

    if config_obj.get(key_name):
        # don't use default paths given in example config file
        if type(config_obj[key_name]) == list and any(path.startswith('path/to/') for path in config_obj[key_name]):
            logging.warning(f'{key_name}: path can not start with "path/to/". Trying to use default instead')
            pass
        elif type(config_obj[key_name]) == str and config_obj[key_name].startswith('path/to/'):
            logging.warning(f'{key_name}: path can not start with "path/to/". Trying to use default instead')
            pass
        else:
            return config_obj[key_name]
    if default: # if default value is given this is evaluated to True
        return default[0]
    else: # key not given in config file and no given default value
        logging.error(f'Error: no default value available for "{key_name}".\
               Please state specifically.')



def pca_cov_variables(vars, include, pbc_indicies):
    """Define coverage variables to be used in PCA.

    Remove coverage variables from the variables used in the PCA if
    coverage_include == FALSE; else append indices for the different
    coverage sets to the cov_vars

    Args:
      vars(str): comma-separated list of variables to be used in PCA
      include(bool): bool if coverage data is included or not
      pbc_indicies: indicies for coverage sets

    Returns:
        (str) modified list of variables to be used in PCA
    """

    cov_vars = "c_cov,c_covsd,g_cov,g_covsd,g_covdev_c"
    out_vars = ""

    if include == 'TRUE':
        for var in vars.split(','):
            if var in cov_vars:
                # add coverage variable with numbered suffix for each
                # coverage data set
                for i in pbc_indicies:
                    out_vars = out_vars + ',' + var + '_' + str(i)
            else:
                out_vars = out_vars + ',' + var
    else:
        for var in vars.split(','):
            if var not in cov_vars:
                out_vars = out_vars + ',' + var

    return out_vars[1:]


def set_config_defaults(config_obj, TAX_DB, db_dir):
    """Set config parameters based on user input and default values.

    Check user input config file and set default values for missing ones;
    write results to temporary config file

    Args:
      config_obj(dict): dict of user input for config parameters

    Returns:
      (dict) dict of processed config parameters with default values
    """

    config_vars = {}

    default_pca_vars = ("c_name,c_num_of_genes,c_len,c_genelenm,c_genelensd," +
                        "g_len,g_lendev_c,g_abspos,g_terminal,c_pearson_r," +
                        "g_pearson_r_o,g_pearson_r_c,c_cov,c_covsd,g_cov," +
                        "g_covsd,g_covdev_c")

    ### check user input config file and set default values for missing ones
    ## General options
    config_vars['fasta_path'] = set_default(config_obj, 'fasta_path')
    config_vars['gff_path'] = set_default(config_obj, 'gff_path')
    config_vars['output_path'] = check_dir_slash(set_default(config_obj,
                                                             'output_path'))
    config_vars['taxon_id'] = int(set_default(config_obj, 'taxon_id'))
    config_vars['taxon_id_rank'] = TAX_DB.taxid2rank[config_vars['taxon_id']]
    config_vars['taxon_name'] = set_default(config_obj, 'taxon_name',
                                            TAX_DB.taxid2name[
                                                config_vars['taxon_id']])

    gff_filename = '.'.join(config_vars['gff_path'].split('/')[-1].split('.')[:-1])
    config_vars['gff_ta_path'] = set_default(config_obj, 'gff_ta_path',
                                             config_vars.get(
                                                 'output_path') + gff_filename + '_w_assigned_taxa.gff')
    config_vars['threads'] = set_default(config_obj, 'threads', 'auto')

    ## Coverage
    # check user input for coverage data and fill missing with default paths
    config_vars['read_paths'] = enumerated_key(config_obj, 'read_paths', [])
    config_vars['bam_paths'] = enumerated_key(config_obj, 'bam_path',
                                              list(config_vars.get('read_paths').keys()),
                                              config_vars.get('output_path')+'mapping_sorted.bam')
    config_vars['pbc_paths'] = enumerated_key(config_obj, 'pbc_path',
                                              list(config_vars.get('bam_paths').keys()),
                                              config_vars.get('output_path')+'pbc.txt')

    # check file existence for each coverage set to determine at which
    # step of the coverage computation the set is
    # dict stores which files to use for next step of computation of PBC
    cov_set_exists = {} # {index: file type to use for futher computation}
    for key in config_vars.get('pbc_paths').keys():
        # reads are required for computation of BAM, BAM for PBC
        # the last file in the list that exists is stored in the dict
        for cov_info in ['read_paths', 'bam_paths', 'pbc_paths']:
            path = config_vars.get(cov_info).get(key)
            if path:
                if check_file_inexistence(path, 'AND') == 'FALSE':
                    # file exists
                    cov_set_exists[key] = cov_info
    config_vars['cov_set_exists'] = cov_set_exists
    # no coverage info available -> no computations
    if len(config_vars.get('cov_set_exists')) == 0:
        config_vars['include_coverage'] = set_default(config_obj,
                                                      'include_coverage',
                                                      'FALSE')
    else:
        config_vars['include_coverage'] = set_default(config_obj,
                                                      'include_coverage',
                                                      'TRUE')
    # if for every coverage set the PBC files exist -> no computations
    if list(config_vars.get('cov_set_exists').values()).count('pbc_paths') == len(config_vars.get('cov_set_exists')) or \
        config_vars.get('include_coverage') == 'FALSE':
        config_vars['compute_coverage'] = set_default(config_obj,
                                                      'compute_coverage',
                                                      'FALSE')
    else:
        config_vars['compute_coverage'] = set_default(config_obj,
                                                      'compute_coverage',
                                                      'TRUE')

    if config_vars.get('include_coverage') == 'FALSE':
        config_vars['pbc_paths'] = {}

    config_vars['min_insert'] = enumerated_key(config_obj,
                                               'min_insert',
                                               list(config_vars.get('read_paths').keys()), '0')
    config_vars['max_insert'] = enumerated_key(config_obj,
                                               'max_insert',
                                               list(config_vars.get('read_paths').keys()), '500')
    config_vars['read_orientation'] = enumerated_key(config_obj,
                                                     'read_orientation',
                                                     list(config_vars.get('read_paths').keys()), 'fr')

    ## Taxonomic assignment
    config_vars['proteins_path'] = set_default(config_obj, 'proteins_path',
                                               config_vars.get(
                                                   'output_path') + 'proteins.faa')
    config_vars['extract_proteins'] = set_default(config_obj,
                                                  'extract_proteins',
                                                  check_file_inexistence(
                                                      config_vars.get(
                                                          'proteins_path'),
                                                      'AND'))
    config_vars['use_phase_info'] = set_default(config_obj, 'use_phase_info',
                                                'auto')
    config_vars['assignment_mode'] = set_default(config_obj, 'assignment_mode',
                                                 'exhaustive')
    # make settings for quick taxonomic assignment mode
    if config_vars.get('assignment_mode') == 'quick':
        config_vars['quick_mode_search_rank'] = set_default(config_obj,
                                                            'quick_mode_search_rank',
                                                            'kingdom')
        config_vars['quick_mode_match_rank'] = set_default(config_obj,
                                                           'quick_mode_match_rank',
                                                           'order')
    else:
        config_vars['quick_mode_search_rank'] = set_default(config_obj,
                                                            'quick_mode_search_rank',
                                                            None)
        config_vars['quick_mode_match_rank'] = set_default(config_obj,
                                                           'quick_mode_match_rank',
                                                           None)
    config_vars['tax_assignment_path'] = set_default(config_obj,
                                                     'tax_assignment_path',
                                                     config_vars.get(
                                                         'output_path') + 'taxonomic_hits.txt')
    if config_vars.get('assignment_mode') == 'quick':
        if type(config_vars.get('tax_assignment_path')) != list:
            if len([config_vars.get('tax_assignment_path')]) != 2:
                path_1 = '.'.join(config_vars.get('tax_assignment_path').split('.')[:-1]) + "_1.txt"
                path_2 = '.'.join(config_vars.get('tax_assignment_path').split('.')[:-1]) + "_2.txt"
                config_vars['tax_assignment_path'] = [path_1, path_2]
    config_vars['compute_tax_assignment'] = set_default(config_obj,
                                                        'compute_tax_assignment',
                                                        check_file_inexistence(
                                                            config_vars.get(
                                                                'tax_assignment_path'),
                                                            'AND'))
    config_vars['database_path'] = set_default(config_obj,
                                               'database_path', f"{db_dir}/nr.dmnd")
    config_vars['taxon_exclude'] = set_default(config_obj, 'taxon_exclude',
                                               'TRUE')
    # use taxon_id_rank for exclusion, if not more specific than "species"
    query_taxon = taxopy.Taxon(config_vars['taxon_id'], TAX_DB)
    if "species" in query_taxon.rank_name_dictionary.keys():
        exclusion_rank = "species"
    else:
        exclusion_rank = config_vars['taxon_id_rank']
    config_vars['exclusion_rank'] = set_default(config_obj, 'exclusion_rank',
                                                exclusion_rank)

    ## Plotting
    config_vars['update_plots'] = set_default(config_obj, 'update_plots',
                                              'FALSE')
    config_vars['num_groups_plot'] = set_default(config_obj, 'num_groups_plot',
                                                 '25')
    config_vars['merging_labels'] = set_default(config_obj, 'merging_labels',
                                                '')
    config_vars['output_pdf'] = set_default(config_obj, 'output_pdf', 'TRUE')
    config_vars['output_png'] = set_default(config_obj, 'output_png', 'FALSE')

    ## Gene info
    config_vars['include_pseudogenes'] = set_default(config_obj,
                                                     'include_pseudogenes',
                                                     'FALSE')
    gff_dict = set_gff_parsing_rules(config_obj)
    config_vars.update(gff_dict)

    ## PCA
    config_vars['input_variables'] = pca_cov_variables(set_default(config_obj,
                                                                   'input_variables',
                                                                   default_pca_vars), config_vars.get('include_coverage'), config_vars.get('pbc_paths').keys())
    config_vars['perform_parallel_analysis'] = set_default(config_obj,
                                                           'perform_parallel_analysis',
                                                           'FALSE')
    config_vars['num_pcs'] = set_default(config_obj, 'num_pcs', '3')
    config_vars['coverage_cutoff_mode'] = set_default(config_obj,
                                                      'coverage_cutoff_mode',
                                                      'default')

    ## Clustering
    config_vars['perform_kmeans'] = set_default(config_obj, 'perform_kmeans',
                                                'FALSE')
    config_vars['kmeans_k'] = set_default(config_obj, 'kmeans_k', 'default')
    config_vars['perform_hclust'] = set_default(config_obj, 'perform_hclust',
                                                'FALSE')
    config_vars['hclust_k'] = set_default(config_obj, 'hclust_k', 'default')
    config_vars['perform_mclust'] = set_default(config_obj, 'perform_mclust',
                                                'FALSE')
    config_vars['mclust_k'] = set_default(config_obj, 'mclust_k', 'default')
    config_vars['perform_dbscan'] = set_default(config_obj, 'perform_dbscan',
                                                'FALSE')
    config_vars['dbscan_groups'] = set_default(config_obj, 'dbscan_groups',
                                               'default')
    config_vars['custom_eps'] = set_default(config_obj, 'custom_eps', '0.3')
    config_vars['custom_minPts'] = set_default(config_obj, 'custom_minPts',
                                               '10')

    return config_vars


def retrieve_tool_paths(script_dir):

    path_dict = eval(open(f"{script_dir}/pathconfig.txt", 'r').readline().strip())
    return_dict = path_dict.get('cmds')
    return_dict['db_dir'] = path_dict.get('data_path')

    return return_dict


def write_cfg2file(config_vars):
    """Write processed config information to yaml file in tmp dir.

    Args:
      config_vars(dict): processed config {parameter name : value}
    """

    with open(config_vars.get('output_path')+'tmp/tmp.cfg.yml', 'w') as out_cfg:
        for key, value in config_vars.items():
                out_cfg.write(f'{key}: {value}\n')


def write_run_overview(config_path, config_vars):
    """Print relevant config information to console.

    Args:
      config_path(str): path to user config file
      config_vars(dict): processed config {parameter name : value}
    """

    logging.info('')
    logging.info(f"Config:\t{config_path}")
    logging.info(f"FASTA:\t{config_vars.get('fasta_path')}")
    logging.info(f"GFF:\t{config_vars.get('gff_path')}")
    logging.info(f"Taxon ID:\t{config_vars.get('taxon_id')}")
    logging.info(f"Output:\t{config_vars.get('output_path')}\n")

    if config_vars.get('include_coverage') == 'TRUE':
        for cov_set, pbc_path in config_vars.get('pbc_paths').items():
            if config_vars.get('cov_set_exists').get(cov_set) == 'pbc_paths':
                logging.info(f"PBC {cov_set}:\t{pbc_path} [exists]")
            else:
                if config_vars.get('cov_set_exists').get(cov_set) == 'bam_paths':
                    basis = 'BAM'
                else:
                    basis = 'read mapping'
                logging.info(f"PBC {cov_set}:\t{pbc_path}")
                logging.info(
                    f"  computation based on:\t"
                    f"{config_vars.get(config_vars.get('cov_set_exists').get(cov_set)).get(cov_set)} "
                    f"[{basis}]")
        logging.info('')
    else:
        logging.info("no valid coverage information provided\n")

    if config_vars.get('extract_proteins') == 'FALSE':
        logging.info(f"Proteins:\t{config_vars.get('proteins_path')} [exists]")
    else:
        logging.info(f"Proteins:\t{config_vars.get('proteins_path')}")
    if config_vars.get('assignment_mode') == 'quick':
        logging.info("Quick assignment mode selected")
        logging.info(f"Filtering search performed "
                     f"on level {config_vars.get('quick_mode_search_rank')}")
        logging.info(f"Hits accepted "
                     f"on level {config_vars.get('quick_mode_match_rank')}")
        if config_vars.get('compute_tax_assignment') == 'FALSE':
            logging.info(f"Taxonomic hits files [exist]:\n"
                         f"{config_vars.get('tax_assignment_path')[0]}"
                         f"{config_vars.get('tax_assignment_path')[1]}")
        else:
            logging.info(f"Taxonomic hits files:\n"
                         f"{config_vars.get('tax_assignment_path')[0]}"
                         f"{config_vars.get('tax_assignment_path')[1]}")

    else:
        logging.info("Exhaustive assignment mode selected")
        if config_vars.get('compute_tax_assignment') == 'FALSE':
            logging.info(f"Taxonomic hits file:"
                         f"{config_vars.get('tax_assignment_path')} [exists]")
        else:
            logging.info(f"Taxonomic hits file:"
                         f"{config_vars.get('tax_assignment_path')}")

    logging.info(f"Query taxon excluded:\t{config_vars.get('taxon_exclude')}\n")

    logging.info(f"Pseudogenes included:\t"
                 f"{config_vars.get('include_pseudogenes')}")
    if config_vars.get('gff_source') != 'default':
        logging.info(f"Rule for GFF parsing:\t{config_vars.get('gff_source')}")
    logging.info('\n')

    logging.info(f"PCA variables:\t{config_vars.get('input_variables')}")
    logging.info(f"Parallel analysis performed:\t"
                 f"{config_vars.get('perform_parallel_analysis')}")
    if config_vars.get('coverage_cutoff_mode') != 'default':
        logging.info(f"Coverage cutoff:\t{config_vars.get('coverage_cutoff_mode')}")

    for clustering in ['kmeans', 'hclust', 'mclust']:
        if config_vars.get('perform_'+clustering) == 'TRUE':
            logging.info(f"{clustering} clustering performed with number of "
                         f"groups:\t{config_vars.get(clustering+'k')}")
    if config_vars.get('perform_dbscan') == 'TRUE':
        logging.info(f"DBSCAN clustering performed with settings:"
                     f"\t{config_vars.get('dbscan_groups')}")
    logging.info('')


def process_config(config_path, script_dir, TAX_DB):
    """Process the users config file and set defaults.

    Processing of user config file, setting defaults, printing overview
    and write processed config parameters to temporary config file

    Args:
      config_path(str): path to user config file
      script_dir(str): path of scripts directory
    """
    # read parameters from config file
    config_obj = yaml.safe_load(open(config_path,'r'))

    cmd_dict = retrieve_tool_paths(script_dir)
    config_vars = set_config_defaults(config_obj, TAX_DB, cmd_dict.get('data_path'))
    config_vars["script_dir"] = script_dir
    config_vars["usr_cfg_path"] = config_path
    config_vars["cfg_path"] = config_vars.get('output_path')+'tmp/tmp.cfg.yml'


    config_vars.update(cmd_dict)

    # write config to file and to console
    write_cfg2file(config_vars)
    write_run_overview(config_path, config_vars)


def make_checks(config_obj):
    """Make data validity checks.

    Cross checks scaffold IDs in FASTA and GFF

    Args:
      config_obj(obj): Config class object of config parameters
    """
    check_assembly_ids(config_obj)
    check_pbc_ids(config_obj)


def cfg2obj(config_path):
    """Make Config class object of config parameters from config file.

    Args:
      config_path(str): path to temporary processed config file

    Returns:
      (obj) Config class object with parameters
    """
    config_obj = yaml.safe_load(open(config_path,'r'))
    cfg = Config(config_obj)

    return cfg


def main():
    """Call module directly with preprocessed config file"""

    config_path = sys.argv[1]
    script_dir = sys.argv[2]

    # read config file and set defaults
    # print to file and summary to console
    process_config(config_path, script_dir)

    # check input data
    make_checks(yaml.safe_load(open(config_path,'r')))


if __name__ == '__main__':
    main()
