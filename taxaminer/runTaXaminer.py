#!/usr/bin/env python

"""taXaminer main script

Runs the pipeline by subsequently calling the submodules. Reads
path to raw config file, extracts and prepares information on this
basis and calls necessary processing steps.

Expects path to raw config file
"""

__author__ = "Freya Arthen"

from . import reduceDims
from . import createOutput
from . import checkInput
from . import compFeatures
from . import compTaxonomicAssignment
from . import parseGFF
from . import createContigOverview

import os
import sys
import yaml
import pathlib
import shutil
import logging
import time
import argparse
import taxopy
import pandas as pd

# custom logging format per level
class Formatter(logging.Formatter):
    """Object for logging."""
    def format(self, record):
        """Set different logging styles for different log levels."""
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        elif record.levelno == logging.WARNING or record.levelno == logging.ERROR:
            self._style._fmt = "%(levelname)s: %(message)s"
        else: #debug
            self._style._fmt = "%(levelname)s [%(asctime)s][%(funcName)s:%(lineno)d]: %(message)s"
        self.datefmt = "%Y-%m-%d %H:%M:%S"
        return super().format(record)


def main():

    # location of script on system to run from outside of directory
    SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

    start_time = time.time()


    # TODO: add more options to pass config parameters via command line
    # TODO: add levels of cleaning (temporarily) created data

    # process logging level information from command line input
    parser = argparse.ArgumentParser()
    parser.add_argument("config_path", help="path to config file")
    parser.add_argument("-q", "--quiet", help="print only warnings and error messages",
                    action="store_true", default=False)
    parser.add_argument("-v", "--verbose", help="print detailed log",
                    action="store_true", default=False)
    parser.add_argument("-k", "--keep", help="keep temporary files",
                    action="store_false", default=True)
    args = parser.parse_args()

    log_level = logging.INFO # default log level
    if args.quiet:
        log_level = logging.WARNING
    elif args.verbose:
        log_level = logging.DEBUG
    # initiate the logger
    logger = logging.getLogger()
    handler = logging.StreamHandler()
    handler.setFormatter(Formatter())
    logger.setLevel(log_level)
    logger.addHandler(handler)

    # reduce verbosity of kaleido/chromium logger when plotting
    logging.getLogger("kaleido").setLevel(logging.ERROR) 
    logging.getLogger("choreographer").setLevel(logging.ERROR) 

    # read config file and retrieve output path
    if pathlib.Path(args.config_path).is_file():
        user_config_dict = yaml.safe_load(open(args.config_path, 'r'))
    else:
        logging.error(f'Config file not found. Please check path:\n\t{args.config_path}')
        sys.exit()
    output_path = user_config_dict.get('output_path')
    if not output_path.endswith('/'): # sanity check path to dir
        output_path = output_path + '/'
    # create dictionaries for output
    pathlib.Path(f'{output_path}tmp/').mkdir(parents=True, exist_ok=True)

    database_dir = eval(
        open(f"{SCRIPT_DIR}/pathconfig.txt", 'r').readline().strip()).get('data_path')

    TAX_DB = taxopy.TaxDb(nodes_dmp=database_dir + "/nodes.dmp",
                          names_dmp=database_dir + "/names.dmp",
                          merged_dmp=database_dir + "/merged.dmp",
                          keep_files=True)

    ## create config file to be used by subsequent modules
    # includes default for parameters not given by user and checks
    # file existence of required files and files to be computed in pipeline
    checkInput.process_config(args.config_path, SCRIPT_DIR, TAX_DB)
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(output_path+'tmp/tmp.cfg.yml')

    pre_time = time.time()
    logging.info('> parsing GFF file')
    gff_df = parseGFF.process(cfg)
    logging.debug(f'finished [{int(time.time() - pre_time)}s]\n')

    if not os.path.isfile(f"{cfg.output_path}gene_info/imputed_gene_table.csv") or cfg.force:

        pathlib.Path(output_path + 'gene_info/').mkdir(parents=True,
                                                       exist_ok=True)

        ## make data integrity checks
        # check GFF and assembly FASTA for ID compatibility
        #checkInput.check_assembly_ids(cfg)

        pre_time = time.time()
        logging.info('> computing gene features')
        compFeatures.process_gene_info(cfg, gff_df)
        logging.debug(f'finished [{int(time.time()-pre_time)}s]\n')

    else:
        logging.info(f"> Reusing gene features in"
                     f" {cfg.output_path}gene_info/imputed_gene_table.csv\n"
                     f"Use option 'force: True' to overwrite.")


    pre_time = time.time()
    logging.info('> executing PCA')
    pathlib.Path(output_path + 'PCA/').mkdir(parents=True, exist_ok=True)
    pca_obj, pca_coordinates, variables = reduceDims.compute_pca(cfg)
    logging.debug(f'finished [{int(time.time()-pre_time)}s]\n')


    # create symlink if proteins.faa file does not exist in taxaminer
    # report (for dashboard)
    if not pathlib.Path(cfg.output_path+'proteins.faa').exists():
        pathlib.Path(cfg.output_path+'proteins.faa').symlink_to(pathlib.Path(cfg.proteins_path).resolve())

    if not os.path.isfile(
            f"{cfg.output_path}taxonomic_assignment/gene_table_taxon_assignment.csv") or cfg.force:
        pre_time = time.time()
        logging.info('> running taxonomic assignment')
        taxonomic_assignment, query_label = compTaxonomicAssignment.run_assignment(cfg, gff_df, pca_coordinates, TAX_DB)
        logging.debug(f'finished [{int(time.time()-pre_time)}s]\n')
    else:
        logging.info(f"> Reusing taxonomic assignment in"
                     f" {cfg.output_path}taxonomic_assignment/gene_table_taxon_assignment.csv\n"
                     f"Use option 'force: True' to overwrite.")
        taxonomic_assignment = pd.read_csv(f"{cfg.output_path}taxonomic_assignment/gene_table_taxon_assignment.csv",
            sep=',', index_col='g_name',
            usecols = ['g_name','fasta_header', 'best_hit', 'best_hitID',
                        'bh_pident', 'bh_evalue', 'bh_bitscore',
                        'lca', 'lcaID', 'refined_lca', 'refined_lcaID',
                        'taxon_assignment', 'taxon_assignmentID', 'is_target',
                        'plot_label', 'plot_labelID'],
            dtype={'taxon_assignmentID': 'str', 'best_hitID': 'str',
                    'refinedLCAID': 'str', 'lcaID': 'str', 'plot_labelID': 'Int64'},
            keep_default_na=False)
        # replace empty value with None
        taxonomic_assignment = taxonomic_assignment.replace(r'^\s*$', None, regex=True)
        target_taxon = compTaxonomicAssignment.init_taxopyon(cfg.taxon_id, set(), TAX_DB)
        query_label = compTaxonomicAssignment.ident_query_label(
                    taxonomic_assignment, target_taxon, TAX_DB)


    pre_time = time.time()
    logging.info('> creating plots')
    all_data_df = createOutput.create_plots(cfg, taxonomic_assignment, pca_obj,
                                            variables, pca_coordinates,
                                            query_label, TAX_DB, gff_df)
    logging.debug(f'finished [{int(time.time() - pre_time)}s]\n')

    pre_time = time.time()
    logging.info('> creating contig summary')
    createContigOverview.process_assignments(cfg, gff_df, all_data_df, TAX_DB)
    logging.debug(f'finished [{int(time.time() - pre_time)}s]\n')


    # make HTML file self-contained, text selectable and change title
    createOutput.perform_adjustments(cfg)

    if args.keep:
        try:
            pre_time = time.time()
            logging.info('> deleting temporary files')
            shutil.rmtree(f"{cfg.output_path}tmp/")
            logging.debug(f'finished [{int(time.time()-pre_time)}s]\n')
        except OSError as e:
            print(f"Error: {cfg.output_path}tmp/ : {e.strerror}")


    logging.info(f'Total runtime: {int(time.time()-start_time)}s\n')


if __name__ == '__main__':
    main()
