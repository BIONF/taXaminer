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
import multiprocessing as mp


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

        self.bam_paths = cfg_dict.get('bam_paths')
        self.include_coverage = cfg_dict.get('include_coverage')

        self.proteins_path = cfg_dict.get('proteins_path')
        self.extract_proteins = cfg_dict.get('extract_proteins')
        self.assignment_mode = cfg_dict.get('assignment_mode')
        self.quick_mode_search_rank = cfg_dict.get('quick_mode_search_rank')
        self.quick_mode_match_rank = cfg_dict.get('quick_mode_match_rank')
        self.tax_assignment_path = cfg_dict.get('tax_assignment_path')
        self.diamond_results_path = cfg_dict.get('diamond_results_path')
        self.slim_diamond_results = cfg_dict.get('slim_diamond_results')
        self.compute_tax_assignment = cfg_dict.get('compute_tax_assignment')
        self.db_dir = f"{cfg_dict.get('db_dir')}/"
        self.database_path = cfg_dict.get('database_path')
        self.target_exclude = cfg_dict.get('target_exclude')
        self.exclusion_rank = cfg_dict.get('exclusion_rank')
        self.ta_generalisation_rank = cfg_dict.get('ta_generalisation_rank')

        self.update_plots = cfg_dict.get('update_plots')
        self.num_groups_plot = cfg_dict.get('num_groups_plot')
        self.merging_labels = cfg_dict.get('merging_labels')
        self.output_pdf = cfg_dict.get('output_pdf')
        self.output_png = cfg_dict.get('output_png')

        self.include_pseudogenes = cfg_dict.get('include_pseudogenes')
        self.input_variables = cfg_dict.get('input_variables')
        self.perform_parallel_analysis = cfg_dict.get('perform_parallel_analysis')
        self.num_pcs = cfg_dict.get('num_pcs')
        self.coverage_cutoff_mode = cfg_dict.get('coverage_cutoff_mode')

        self.script_dir = cfg_dict.get('script_dir')
        self.usr_cfg_path = cfg_dict.get('usr_cfg_path')
        self.cfg_path = cfg_dict.get('cfg_path')

        # tool cmds
        self.diamond = cfg_dict.get('diamond')
        self.samtools = cfg_dict.get('samtools')
        self.bedtools = cfg_dict.get('bedtools')
        self.bowtie2 = cfg_dict.get('bowtie2')
        self.krona = cfg_dict.get('krona')


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


def enumerated_key(config_obj, key_name):
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
    keys = [x for x in range(0,50)]
    for match in matches:
        mod_default = None
        if config_obj.get(match):
            if len(match.split('_')) >= 3 and match.split('_')[2].isdigit():
                match_num = int(match.split('_')[2])
                dict[match_num] = set_default(config_obj, match)
                keys.remove(match_num)
            else:
                dict[keys[0]] = set_default(config_obj, match)
                del keys[0]

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
        return False
    else:
        return True


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
        return config_obj[key_name]
    if default: # if default value is given this is evaluated to True
        return default[0]
    else: # key not given in config file and no given default value
        logging.error(f'Error: no default value available for "{key_name}".\
               Please state specifically.')


def set_yesno_default(config_obj, key_name, *default):
    """Check if parameter is given in config file, else set to default.

    Args:
      config_obj(dict): dict of user input for config parameters
      key_name(str): key to look for in config_obj
      *default(str): default parameter for key_name

    Returns:
        (str) if key_name in config_obj then value, else default

    """

    if config_obj.get(key_name):
        if config_obj.get(key_name).lower() == "true" or \
            config_obj.get(key_name).lower() == "yes" or \
            config_obj.get(key_name).lower() == "y":
            return True
        else:
            return False
    if default: # if default value is given this is evaluated to True
        return default[0]
    else: # key not given in config file and no given default value
        logging.error(f'Error: no default value available for "{key_name}".\
               Please state specifically.')



def pca_cov_variables(vars, include, bam_indicies):
    """Define coverage variables to be used in PCA.

    Remove coverage variables from the variables used in the PCA if
    coverage_include == FALSE; else append indices for the different
    coverage sets to the cov_vars

    Args:
      vars(str): comma-separated list of variables to be used in PCA
      include(bool): bool if coverage data is included or not
      bam_indicies: indicies for coverage sets

    Returns:
        (str) modified list of variables to be used in PCA
    """

    cov_vars = "c_cov,c_covsd,g_cov,g_covsd,g_covdev_c"
    out_vars = ""

    if include:
        for var in vars.split(','):
            if var in cov_vars:
                # add coverage variable with numbered suffix for each
                # coverage data set
                for i in bam_indicies:
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
    config_vars['threads'] = set_default(config_obj, 'threads', mp.cpu_count())

    ## Coverage
    config_vars['bam_paths'] = enumerated_key(config_obj, 'bam_path')
    config_vars['include_coverage'] = set_yesno_default(config_obj, 'include_coverage',
                                                  True if config_vars['bam_paths'] else False)

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
    config_vars['diamond_results_path'] = set_default(config_obj,
                                                     'diamond_results_path',
                                                     config_vars.get(
                                                         'output_path') + 'tmp/diamond_results.txt')
    config_vars['slim_diamond_results'] = set_yesno_default(config_obj,
                                                       'slim_diamond_results',
                                                       False)

    if config_vars.get('assignment_mode') == 'quick':
        if type(config_vars.get('tax_assignment_path')) != list:
            if len([config_vars.get('tax_assignment_path')]) != 2:
                path_1 = '.'.join(config_vars.get('tax_assignment_path').split('.')[:-1]) + "_1.txt"
                path_2 = '.'.join(config_vars.get('tax_assignment_path').split('.')[:-1]) + "_2.txt"
                config_vars['tax_assignment_path'] = [path_1, path_2]
    config_vars['compute_tax_assignment'] = set_yesno_default(config_obj,
                                                        'compute_tax_assignment',
                                                        check_file_inexistence(
                                                            config_vars.get(
                                                                'tax_assignment_path'),
                                                            'AND'))
    config_vars['database_path'] = set_default(config_obj,
                                               'database_path', f"{db_dir}/nr.dmnd")
    config_vars['target_exclude'] = set_yesno_default(config_obj, 'target_exclude',
                                               True)
    # use taxon_id_rank for exclusion, if not more specific than "species"
    query_taxon = taxopy.Taxon(config_vars['taxon_id'], TAX_DB)
    if "species" in query_taxon.rank_name_dictionary.keys():
        exclusion_rank = "species"
    else:
        exclusion_rank = config_vars['taxon_id_rank']
    config_vars['exclusion_rank'] = set_default(config_obj, 'exclusion_rank',
                                                exclusion_rank)
    config_vars['ta_generalisation_rank'] = set_default(config_obj, 'ta_generalisation_rank',
                                                'genus')

    ## Plotting
    config_vars['update_plots'] = set_yesno_default(config_obj, 'update_plots',
                                              False)
    config_vars['num_groups_plot'] = set_default(config_obj, 'num_groups_plot',
                                                 '25')
    config_vars['merging_labels'] = set_default(config_obj, 'merging_labels',
                                                '')
    config_vars['output_pdf'] = set_yesno_default(config_obj, 'output_pdf', True)
    config_vars['output_png'] = set_yesno_default(config_obj, 'output_png', False)

    ## Gene info
    config_vars['include_pseudogenes'] = set_yesno_default(config_obj,
                                                     'include_pseudogenes',
                                                     False)

    ## PCA
    config_vars['input_variables'] = pca_cov_variables(set_default(config_obj,
                                                                   'input_variables',
                                                                   default_pca_vars),
                                                       config_vars.get('include_coverage'),
                                                       config_vars.get('bam_paths').keys())
    config_vars['perform_parallel_analysis'] = set_yesno_default(config_obj,
                                                           'perform_parallel_analysis',
                                                           False)
    config_vars['num_pcs'] = set_default(config_obj, 'num_pcs', '3')
    config_vars['coverage_cutoff_mode'] = set_default(config_obj,
                                                      'coverage_cutoff_mode',
                                                      'default')

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

    if config_vars.get('include_coverage'):
        for cov_set, bam_path in config_vars.get('bam_paths').items():
            logging.info(f"BAM file {cov_set}:\t{bam_path}\n")
    else:
        logging.info("no valid coverage information provided\n")

    if not config_vars.get('extract_proteins'):
        logging.info(f"Proteins:\t{config_vars.get('proteins_path')} [exists]")
    else:
        logging.info(f"Proteins:\t{config_vars.get('proteins_path')}")
    if config_vars.get('assignment_mode') == 'quick':
        logging.info("Quick assignment mode selected")
        logging.info(f"Filtering search performed "
                     f"on level {config_vars.get('quick_mode_search_rank')}")
        logging.info(f"Hits accepted "
                     f"on level {config_vars.get('quick_mode_match_rank')}")
        if not config_vars.get('compute_tax_assignment'):
            logging.info(f"Taxonomic hits files [exist]:\n"
                         f"{config_vars.get('tax_assignment_path')[0]}"
                         f"{config_vars.get('tax_assignment_path')[1]}")
        else:
            logging.info(f"Taxonomic hits files:\n"
                         f"{config_vars.get('tax_assignment_path')[0]}"
                         f"{config_vars.get('tax_assignment_path')[1]}")

    else:
        logging.info("Exhaustive assignment mode selected")
        if not config_vars.get('compute_tax_assignment'):
            logging.info(f"Taxonomic hits file:"
                         f"{config_vars.get('tax_assignment_path')} [exists]")
        else:
            logging.info(f"Taxonomic hits file:"
                         f"{config_vars.get('tax_assignment_path')}")

    logging.info(f"Query taxon excluded:\t{config_vars.get('target_exclude')}\n")

    #TODO: are pseudogenes still supported?
    logging.info(f"Pseudogenes included:\t"
                 f"{config_vars.get('include_pseudogenes')}")
    logging.info('\n')

    logging.info(f"PCA variables:\t{config_vars.get('input_variables')}")
    logging.info(f"Parallel analysis performed:\t"
                 f"{config_vars.get('perform_parallel_analysis')}")
    if config_vars.get('coverage_cutoff_mode') != 'default':
        logging.info(f"Coverage cutoff:\t{config_vars.get('coverage_cutoff_mode')}")

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
    config_vars = set_config_defaults(config_obj, TAX_DB, cmd_dict.get('db_dir'))
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
    #TODO: check_gff_format()


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
