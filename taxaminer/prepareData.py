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

import sys
from Bio import Seq

from . import checkInput


def set_seqs(cfg, gff_df, gene, contig_seq, proteins_file):
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
        return (None, None)

    seq = ''

    # sort coding feature by start coordinate
    sorted_coding_features = sorted(gene.coding_features, key=lambda cds: gff_df.loc[cds].start)
    transl_tables = set()
    for coding_feature in sorted_coding_features:
        cds = gff_df.loc[coding_feature]
        transl_tables.add(str(cds.transl_table))
        if not cfg.use_phase:
            # per default do not use the phase info
            seq += str(Seq.Seq(contig_seq[cds.start-1:cds.end]))
        else:
            if gene.strand == '-':
                seq += str(Seq.Seq(
                    contig_seq[(cds.start - 1):(cds.end- (int(cds.phase) if cds.phase != '.' else 0))]))
            else:
                seq += str(Seq.Seq(contig_seq[(cds.start - 1 + (int(cds.phase) if cds.phase != '.' else 0)):(cds.end)]))

    if gene.strand == '-':
        seq = str(Seq.Seq(seq).reverse_complement())

    if len(seq)%3 != 0:
        seq += ("N"*(3-(len(seq)%3)))

    if len(transl_tables) == 1:
        transl_table = transl_tables.pop()
    if (len(transl_tables) >= 2) or (
            transl_table != gene.transl_table and gene.transl_table != '1'):
        # a translational table is in most cases only given for CDS features
        # gene feature has default value of 1 given for translational table
        logging.warning(f'Ambiguous translational table for gene {gene.name}. '
                        f'(Values given in GFF: {",".join(transl_tables)}{gene.transl_table})'
                        f'Continuing with table given for gene entry.')
        transl_table = gene.transl_table
    protein_seq = str(Seq.Seq(seq).translate(table=int(transl_table)))

    #>#write_seq2file(proteins_file, gene.transcript_id, protein_seq)

    return (gene.transcript_id, protein_seq)


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
                        genes.apply(lambda row: set_seqs(cfg, gff_df, row, contig_seq, proteins_file), axis=1)
                    # retrieve ID from current contig
                    current_contig = line.strip().lstrip(">").split()[0]
                    contig_seq = ''
                else:
                    contig_seq += line.strip()
            # add genes from last contig
            if not gff_df.loc[(gff_df['scaffold'] == current_contig) & (gff_df['type'] == "gene")].empty:
                genes = gff_df.loc[(gff_df['scaffold'] == current_contig) & (gff_df['type'] == "gene")]
                genes.apply(lambda row: set_seqs(cfg, gff_df, row, contig_seq, proteins_file), axis=1)


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
