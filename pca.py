#!/usr/bin/env python

"""compute PCA of gene desprictors



Expects path to config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

from classes import Assembly, Contig, Gene
import prepare_and_check

# import numpy as np
# import pathlib # to create directories
# import operator # for quick comparisons
# import scipy.stats as stats # for Pearson's R
# from itertools import product as itertools_product # to generate all possible oligonucleotides from base alphabet
# from Bio.Seq import Seq as BioPython_Seq # to count oligonucleotides (also overlapping ones! not implemented in normal count)
# import sys # parse command line arguments

import logging
import pandas as pd


def compute_pca(cfg):

    df = pd.read_csv(cfg.output_path+'gene_info/imputed_gene_table.txt', header=0)
    print(df)




def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(config_path)

    compute_pca(cfg)


if __name__ == '__main__':
    main()
