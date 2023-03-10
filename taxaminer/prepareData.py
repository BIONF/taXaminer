#!/usr/bin/env python

"""Prepare data for running taXaminer

#TODO: parse gff file here
Parses the GFF file

Prepares coverage data for further processing. Raw reads are mappend
and mapping files (BAM) are converted to per base coverage (PBC) files

Compute protein FASTA by parsing gene sequences baed on coordinates from
GFF and FASTA with scaffold sequences (assembly FASTA) and converting
it to amino acid sequence

Expects prepared config file with information regarding available
coverage information and corresponding paths to find the files
(preparation of config by prepare_and_check.py)
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

import logging
import sys
import pathlib
import subprocess
from Bio import Seq

from . import checkInput


def set_seqs(gff_df, gene, contig_seq, proteins_file):
    """Extract nuc sequences and translate to AA for genes on contig.

    Use coordinates of longest transcript in gene Feature object
    to extract the gene sequence and translate it into AA sequence

    Args:
      proteins_file(obj): file object for protein FASTA file
      contigs(dict): {scaffold ID : [gene IDs]}
      features(dict): {GFF feature ID : Feature object}
      current_contig(str): ID of currrent contig
      contig_seq(str): sequence for current_contig
    """

    if not gene.coding_features:
        # no transcript was found for gene
        return

    seq = ''
    for coding_feature in gene.coding_features:
        cds = gff_df.loc[coding_feature]
        # on reverse strand CDS needs to be reverse complemented before concatenated
        if gene.strand == '-':
            seq += str(Seq.Seq(contig_seq[cds.start-1:cds.end]).reverse_complement())
        else:
            seq += str(Seq.Seq(contig_seq[cds.start-1:cds.end]))


    if len(seq)%3 != 0:
        seq += ("N"*(3-(len(seq)%3)))

    protein = str(Seq.Seq(seq).translate(table=int(gene.transl_table)))

    write_seq2file(proteins_file, gene.transcript_id, protein)


def write_seq2file(file, seq_id, seq):
    """Write proteins sequences to file

    Args:
      proteins_file(obj): file object for protein FASTA file
      contigs(dict): {scaffold ID : [gene IDs]}
      features(dict): {GFF feature ID : Feature object}
      current_contig(str): ID of currrent contig
      contig_seq(str): sequence for current_contig
    """

    file.write(">" + seq_id + "\n")
    for i in range(0,len(seq),80):
        file.write(seq[i:i+80]+"\n")


def extract_seq(cfg, gff_df):
    """Retrieve sequence for contig from FASTA and write protein FASTA.

    Read assembly FASTA, retrieve scaffold sequence and call set_seqs()
    to write AA sequences for genes to file

    Args:
      cfg(obj): Config class object of config parameters
      contigs(dict): {scaffold ID : [gene IDs]}
      gff_df(dict): {GFF feature ID : Feature object}
    """

    current_contig = ""

    with open(cfg.fasta_path, "r") as fasta_file:
        with open(cfg.proteins_path, "w") as proteins_file:
            for line in fasta_file:
                if line.startswith(">"):
                    # False if no genes in contig dict (geneless contigs)
                    if not gff_df.loc[(gff_df['scaffold'] == current_contig) & (gff_df['type'] == "gene")].empty:
                        genes = gff_df.loc[(gff_df['scaffold'] == current_contig) & (gff_df['type'] == "gene")]
                        genes.apply(lambda row: set_seqs(gff_df, row, contig_seq, proteins_file), axis=1)
                    # retrieve ID from current contig
                    current_contig = line.strip().lstrip(">").split()[0]
                    contig_seq = ''
                else:
                    contig_seq += line.strip()
            # add genes from last contig
            if not gff_df.loc[(gff_df['scaffold'] == current_contig) & (gff_df['type'] == "gene")].empty:
                genes = gff_df.loc[(gff_df['scaffold'] == current_contig) & (gff_df['type'] == "gene")]
                genes.apply(lambda row: set_seqs(gff_df, row, contig_seq, proteins_file), axis=1)


def generate_fasta(cfg, gff_df):
    """Generate FASTA file of protein sequences.

    Extract coordinates for genes and transcripts and their relation,
    determine the longest transcript for each gene,
    extract the nucleotide sequence and translate to amino acid sequence

    Args:
      cfg(obj): Config class object of config parameters
      :param gff_df:
    """

    extract_seq(cfg, gff_df)


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    #TODO: write gff_df to file/pickle
    gff_df = sys.argv[2]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    generate_fasta(cfg, gff_df)


if __name__ == '__main__':
    main()
