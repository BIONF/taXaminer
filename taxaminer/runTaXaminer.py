#!/usr/bin/env python

"""taXaminer main script

Runs the pipeline by subsequently calling the submodules. Reads
path to raw config file, extracts and prepares information on this
basis and calls necessary processing steps.

Expects path to raw config file
"""

__author__ = "Freya Arthen"
__version__ = "0.6.0"

# package specific modules
import pickle

from . import reduceDims
from . import createOutput
from . import checkInput
from . import prepareData
from . import compFeatures
from . import compTaxonomicAssignment
from . import parseGFF
from . import createContigOverview

import os
import sys
import yaml
import pathlib
import subprocess
import shutil
import logging
import time
import argparse
import taxopy


# custom logging format per level
class Formatter(logging.Formatter):
    """Object for logging."""
    def format(self, record):
        """Set different logging styles for different log levels."""
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        else:
            self._style._fmt = "[%(asctime)s] %(levelname)s [%(funcName)s:%(lineno)d] %(message)s"
        self.datefmt = "%Y-%m-%d %H:%M:%S"
        return super().format(record)


def main():
    """Run the program and call submodules."""

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
                    action="store_false", default=False)
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

    # read config file and retrieve output path
    if pathlib.Path(args.config_path).is_file():
        user_config_dict = yaml.safe_load(open(args.config_path, 'r'))
    else:
        logging.error(f'config file not found:   {args.config_path}')
        sys.exit()
    output_path = user_config_dict.get('output_path') # location for results
    if not output_path.endswith('/'): # sanity check path to dir
        output_path = output_path + '/'
    # create dictionaries for output
    pathlib.Path(output_path+'tmp/').mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_path+'gene_info/').mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_path+'PCA/').mkdir(parents=True, exist_ok=True)

    database_dir = eval(open(f"{SCRIPT_DIR}/pathconfig.txt",
                             'r').readline().strip()).get('data_path')

    # TODO: precompute the tax_db during setup
    TAX_DB = taxopy.TaxDb(nodes_dmp=database_dir + "/nodes.dmp",
                          names_dmp=database_dir + "/names.dmp",
                          merged_dmp=database_dir + "/merged.dmp",
                          keep_files=True)

    # f = open('tax_db.obj', 'wb')
    # pickle.dump(TAX_DB, f)
    # f.close()

    # f = open('tax_db.obj', 'rb')
    # TAX_DB = pickle.load(f)
    # f.close()


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

    if not cfg.update_plots:

        ## make data integrity checks
        # check GFF and assembly FASTA for ID compatibility
        #checkInput.check_assembly_ids(cfg)

        pre_time = time.time()
        logging.info('> computing gene descriptors')
        compFeatures.process_gene_info(cfg, gff_df)
        logging.debug(f'finished [{int(time.time()-pre_time)}s]\n')

        pre_time = time.time()
        logging.info('> executing PCA and clustering')
        pca_obj, pca_coordinates, variables = reduceDims.compute_pca(cfg)
        logging.debug(f'finished [{int(time.time()-pre_time)}s]\n')

        # create symlink if proteins.faa file does not exist in taxaminer
        # report (for dashboard)
        if not pathlib.Path(cfg.output_path+'proteins.faa').exists():
            pathlib.Path(cfg.output_path+'proteins.faa').symlink_to(pathlib.Path(cfg.proteins_path).resolve())

    pre_time = time.time()
    logging.info('> running taxonomic assignment')
    taxonomic_assignment, query_label = compTaxonomicAssignment.run_assignment(cfg, gff_df, pca_coordinates, TAX_DB)
    logging.debug(f'finished [{int(time.time()-pre_time)}s]\n')

    pre_time = time.time()
    logging.info('> creating plots')
    all_data_df = createOutput.create_plots(cfg, taxonomic_assignment, pca_obj,
                                            variables, pca_coordinates,
                                            query_label, TAX_DB, gff_df)
    logging.debug(f'finished [{int(time.time() - pre_time)}s]\n')

    # pre_time = time.time()
    # logging.info('> creating contig summary')
    # createContigOverview.process_assignments(cfg, gff_df, all_data_df, TAX_DB)
    # logging.debug(f'finished [{int(time.time() - pre_time)}s]\n')


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
