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

from . import checkData


class Feature:
    """Object for GFF feature
    """
    def __init__(self, contig, info_dict):

        self.contig = contig
        self.id = info_dict.get('ID')
        self.parent = info_dict.get('Parent')
        self.start = info_dict.get('start')
        self.end = info_dict.get('end')
        self.strand = info_dict.get('strand')
        self.phase = info_dict.get('phase')
        self.biotype = info_dict.get('biotype') if info_dict.get('biotype') else info_dict.get('gene_biotype')
        self.transl_table = info_dict.get('transl_table')

        self.children = {}
        self.parsed_parent = None

        # gene only
        self.transcripts = {} # key = transcript ID, value = CDS length
        self.cdss = {} # key = transcript ID, value = [CDS IDs]
        self.coordinates = [] # coordinates of CDS of currently longest transcript
        self.phased_coordinates = []    # coordinates of CDS of currently longest transcript
                                        # adjusted for the information of phase attribute
        self.prot_seqs = {'phased':'', 'non-phased': ''} # store both seqs

def bam2pbc(cfg, bam, pbc):
    """Convert BAM file to PBC.

    Uses bedtools genomecov with option -d to convert mapping
    information in BAM format to per base coverage file (PBC)
    Option -d reports read depth at each position

    Args:
      bam(str): path to BAM file
      pbc(str): path to PBC file (output)
    """

    # FIXME: check this here (add specific output filtering)
    # --> probably error because samtools view command was skipped in mapping()
    cmd = f'{cfg.bedtools} genomecov -ibam "{bam}" -d > "{pbc}"'
    out = subprocess.run([cmd], shell=True, capture_output=True)
    if out.returncode != 0:
        logging.error(f'converting BAM to PBC failed:\n{cmd}')
        logging.error('Error message:\n'+out.stderr.decode())
        sys.exit()


def mapping(cfg, mapping_dir, fasta_path, read_paths, min_insert_set,
            max_insert_set, read_orientation_set, bam_path, threads):
    """Mapping reads to assembly using Bowtie2.

    Uses Bowtie2 to map reads to reference assembly. Options used:
      -a: every alignment is reported (multiple alignments enabled)
      --sensitive: sensitivity mode
      -I/-X: insert size (lower and upper limit)

    Aligns paired and single end reads sets based on how many paths
    there are in read_paths

    Args:
      mapping_dir(str): path to output directory of mapping files
      fasta_path(str): path to FASTA file of assembly
      read_paths(str,list): (list of) path(s) to read file(s)
      insert_size(int): insert size for paired end reads, else None
      bam_path(str): path to BAM file
    """

    # TODO: check if SAM file exists at default location
    if type(read_paths) != list:
        read_paths = [read_paths]

    ## prepare commands
    cmd_build = f'{cfg.bowtie2}-build "{fasta_path}" {mapping_dir}assembly'
    # depending on number of paths in list read_paths different cmd
    # for mapping of paired end or single end data is used
    if len(read_paths) == 2:
        cmd_mapping =  f'{cfg.bowtie2} --sensitive ' \
                       f'-I {min_insert_set} -X {max_insert_set} ' \
                       f'--{read_orientation_set} -a -p {threads} ' \
                       f'-x {mapping_dir}assembly ' \
                       f'-1 {read_paths[0]} ' \
                       f'-2 {read_paths[1]} -S {mapping_dir}mapping.sam'
    elif len(read_paths) == 1:
        cmd_mapping =  f'{cfg.bowtie2} --sensitive -a -p {threads} ' \
                       f'-x {mapping_dir}assembly -U {read_paths[0]} ' \
                       f'-S {mapping_dir}mapping.sam'
    else:
        print(f'Error: unexpected number of read paths detected ' \
              f'for coverage set (max allowed: 2; given: {len(read_paths)}). ' \
              f'Please recheck your input.')
    cmd_view = f'{cfg.samtools} view -b {mapping_dir}mapping.sam ' \
               f'| {cfg.samtools} sort -o {bam_path}'

    ## run commands
    # build index
    out = subprocess.run([cmd_build], shell=True, capture_output=True, check=True)
    if out.returncode != 0:
        logging.error(f'failed building index for mapping:\n{cmd_build}')
        logging.error('Error message:\n'+out.stderr.decode())
        sys.exit()
    # mapping
    out = subprocess.run([cmd_mapping], shell=True, capture_output=True, check=True)
    if out.returncode != 0:
        logging.error(f'failed aligning reads to assembly:\n{cmd_mapping}')
        logging.error('Error message:\n'+out.stderr.decode())
        sys.exit()
    # convert SAM to BAM
    out = subprocess.run([cmd_view], shell=True, capture_output=True, check=True)
    if out.returncode != 0:
        logging.error(f'sorting and converting SAM to BAM failed:\n{cmd_view}')
        logging.error(f'Error message:\n{out.stderr.decode()}')
        sys.exit()

def process_coverage(cfg):
    """Prepare coverage information for further analysis

    Uses information in config object to perform computations for
    further use of coverage information. Calls
      mapping()
      bam2pbc()

    Args:
      cfg: Config object with config parameters
    """

    if cfg.compute_coverage:
        # cfg.cov_set_exists: {index: status}
        # for every set of coverage data there is
        for cov_set, status in cfg.cov_set_exists.items():
            # check the status
            if status == 'pbc_paths':
                # coverage data already in PBC format
                continue

            elif status == 'bam_paths':
                # coverage data in BAM format -> convert to PBC
                bam_path = cfg.bam_paths.get(cov_set)
                pbc_path = cfg.pbc_paths.get(cov_set)
                bam2pbc(cfg, bam_path, pbc_path)

            elif status == 'read_paths':
                # coverage data as read files -> map + convert
                # create directory to store mapping files
                mapping_dir = f'{cfg.output_path}mapping_files_{cov_set}/'
                pathlib.Path(mapping_dir).mkdir(parents=True, exist_ok=True)

                read_paths_set = cfg.read_paths.get(cov_set)
                min_insert_set = int(cfg.min_insert.get(cov_set))
                max_insert_set = cfg.max_insert.get(cov_set)
                read_orientation_set = cfg.read_orientation.get(cov_set)
                bam_path = cfg.bam_paths.get(cov_set)
                pbc_path = cfg.pbc_paths.get(cov_set)

                num_threads = cfg.threads if cfg.threads != 'auto' else 1


                mapping(cfg, mapping_dir, cfg.fasta_path, read_paths_set,
                        min_insert_set, max_insert_set, read_orientation_set, bam_path, num_threads)
                bam2pbc(cfg, bam_path, pbc_path)

            else:
                print(f"Error: status of coverage set {cov_set} unidentified")


def get_cds_coordinates(cfg):
    """ Parse GFF to find gene and their corresponding mRNA and/or CDS

    Parses the GFF and, while considering the special parsing options of
    taXamine, gathers information on genes, mRNA and CDS. Saves for each
    contig the list of genes IDs that are annotated on it. CDS and transcripts
    features are connected to genes. Transcript Feature object are stored
    together with genes in features dictionary; CDS features are saved in
    gene object.

    Can not handle unsorted GFFs!

    Args:
      cfg(obj): Config class object of config parameters

    Returns:
        contigs(dict): {scaffold ID : [gene IDs]}
        features(dict): {GFF feature ID : Feature object}
    """

    contigs = {} #key = contig, value = [gene ID]
    features = {} #key = feature ID, value = feature object

    unmatched_lines = {} # key = parent_id, value = line

    with open(cfg.gff_path, "r") as gff_file:
        parent = None
        for line in gff_file:
            if not line.startswith("#"):
                spline = line.strip().split("\t")
                # only read lines that should be considered based on cfg.gff_source
                # if default is selected, read all lines
                if cfg.gff_source == "default"  or spline[1] == cfg.gff_source:
                    if spline[2] in [cfg.gff_gene_type, cfg.gff_transcript_type,
                                    cfg.gff_cds_type]:
                        #gather information for GFF entries of type gene, mRNA and CDS
                        contig = spline[0]
                        info_dict={"start": int(spline[3]), "end": int(spline[4]), "strand": spline[6], "phase": spline[7]}
                        # save each attribute separately in dict
                        for elem in spline[8].split(";"):
                            # check if attributes are empty/faulty
                            if len(elem.split("=")) == 2:
                                key, value = elem.split("=")
                                info_dict[key] = value
                        feature = Feature(contig, info_dict)
                        if spline[2] == cfg.gff_gene_type: # gene
                            # if biotype tag is available, use only genes with protein_coding tag
                            if (feature.biotype and feature.biotype == "protein_coding") or not feature.biotype:
                                features[feature.id] = feature
                                if contig in contigs.keys():
                                    contigs[contig].append(feature.id)
                                else:
                                    contigs[contig] = [feature.id]
                        else: # mRNA or CDS
                            # get Feature object of parent for current feature
                            if cfg.gff_connection == 'parent/child' or cfg.gff_source in ['default', 'maker', 'augustus_masked']:
                                parent_attr = cfg.gff_parent_attr.lower() if cfg.gff_parent_attr else 'parent'
                                parent_id = getattr(feature,parent_attr)
                                parent = features.get(parent_id)
                            elif cfg.gff_connection == 'inline':
                                parent_id = getattr(feature,cfg.gff_gene_attr.lower())
                                parent = features.get(parent_id)
                            else:
                                logging.error('Please check GFF parsing rule. Unclear how to connect gene and transcript')
                            if parent:
                                # which type is parent feature: gene or transcript?
                                feature.parsed_parent = parent.id
                                if spline[2] == cfg.gff_transcript_type:
                                    # current feature -> transcript
                                    # parent -> gene
                                    features[feature.id] = feature
                                    parent.transcripts[feature.id] = 0
                                    parent.cdss[feature.id] = []
                                else:
                                    # current feature -> CDS
                                    # parent -> gene or transcript
                                    if parent.parsed_parent:
                                        # parsed_parent of parent is already assigned
                                        # parent -> transcript
                                        gene = features.get(parent.parsed_parent)
                                        transcript_id = parent.id
                                    else:
                                        # parent -> gene
                                        # thus transcripts are not specified and
                                        # no transcripts are available
                                        # current feature -> transcript and CDS
                                        gene = parent
                                        transcript_id = gene.id
                                        parent.transcripts[transcript_id] = 0
                                        if not parent.cdss.get(transcript_id):
                                            parent.cdss[transcript_id] = []
                                    gene.cdss[transcript_id].append(feature)
                            else:
                                if parent_id in unmatched_lines.keys():
                                    unmatched_lines[parent_id].append(spline)
                                else:
                                    unmatched_lines[parent_id] = [spline]
            elif "#FASTA" in line: # if FASTA block has been reached
                break

    # GFF3 was unsorted
    count_it = 0
    while unmatched_lines and count_it <= 2:
        # make at max two iterations:
        # genes are already added to features because they need no parent
        # mRNA features can be added in first iteration as their parent (gene)
        # will definetly be there
        # CDS features will be at last added in second one (if not first)
        count_it += 1
        unmatched_lines_copy = {key: value for key, value in unmatched_lines.items()}
        for parent_id, lines in unmatched_lines_copy.items():
            if parent_id in features.keys():
                del unmatched_lines[parent_id]
                for spline in lines:
                    #TODO: move this redundant code block to a function
                    #gather information for GFF entries of type gene, mRNA and CDS
                    contig = spline[0]
                    info_dict={"start": int(spline[3]), "end": int(spline[4]), "strand": spline[6], "phase": spline[7]}
                    # save each attribute separately in dict
                    for elem in spline[8].split(";"):
                        # check if attributes are empty/faulty
                        if len(elem.split("=")) == 2:
                            key, value = elem.split("=")
                            info_dict[key] = value
                    feature = Feature(contig, info_dict)
                    parent = features.get(parent_id)
                    # which type is parent feature: gene or transcript?
                    feature.parsed_parent = parent.id
                    if spline[2] == cfg.gff_transcript_type:
                        # current feature -> transcript
                        # parent -> gene
                        features[feature.id] = feature
                        parent.transcripts[feature.id] = 0
                        parent.cdss[feature.id] = []
                    else:
                        # current feature -> CDS
                        # parent -> gene or transcript
                        if parent.parsed_parent:
                            # parsed_parent of parent is already assigned
                            # parent -> transcript
                            gene = features.get(parent.parsed_parent)
                            transcript_id = parent.id
                        else:
                            # parent -> gene
                            # thus transcripts are not specified and
                            # no transcripts are available
                            # current feature -> transcript and CDS
                            gene = parent
                            transcript_id = gene.id
                            parent.transcripts[transcript_id] = 0
                            if not parent.cdss.get(transcript_id):
                                parent.cdss[transcript_id] = []
                        gene.cdss[transcript_id].append(feature)

    return contigs, features

def get_longest_transcript(contigs, features):
    """Identify the transcript with the longest CDS for each gene.

    Computes length of each transcript associated with gene Features
    and stores the coordinates for the respective CDS. Also corrects
    for phase, adds strand and translational table information

    Args:
      contigs(dict): {scaffold ID : [gene IDs]}
      features(dict): {GFF feature ID : Feature object}
    """

    for contig, c_genes in contigs.items():
        cds = None
        for gene_id in c_genes:
            gene = features.get(gene_id)

            if not gene.transcripts or not gene.cdss:
                # no transcript was found for gene
                continue

            for transcript_id, cdss in gene.cdss.items():
                length = 0
                for cds in cdss:
                    length += cds.end-cds.start+1
                gene.transcripts[transcript_id] = length
            # identify transcript with longest CDS
            max_len_t = max(gene.transcripts, key=gene.transcripts.get)
            # get coordinates for that CDS
            for cds in gene.cdss.get(max_len_t):
                gene.coordinates.append((cds.start-1,cds.end,cds.phase))

            if not cds or not gene.coordinates:
                continue

            if cds.strand == '+':
                gene.strand = False # CDS is on forward strand
            else:
                gene.strand = True


            # gene.coordinates[0] = (gene.coordinates[0][0]+(int(gene.coordinates[0][2]) if gene.coordinates[0][2] != '.' else 0),
            #        gene.coordinates[0][1], gene.coordinates[0][2])

            # correct coordinates for phase
            if gene.strand: # minus strand
                gene.phased_coordinates = [(coordinate[0],
                                (coordinate[1]-(int(coordinate[2]) if coordinate[2] != '.' else 0)), coordinate[2]) for coordinate in gene.coordinates]
            if not gene.strand: # plus strand
                gene.phased_coordinates = [((coordinate[0]+(int(coordinate[2]) if coordinate[2] != '.' else 0)),
                                coordinate[1], coordinate[2]) for coordinate in gene.coordinates]
            # sort coordinates
            gene.coordinates = sorted(gene.coordinates, key=lambda k: k[0], reverse=gene.strand)


            # add translational table information if available
            gene.transl_table = cds.transl_table
            if not gene.transl_table:
                gene.transl_table = '1'

            # print(gene.__dict__)


def set_seqs(contigs, features, current_contig, contig_seq,stop_codons, outta_frame):
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



    for gene_id in contigs.get(current_contig):
        gene = features.get(gene_id)
        # print(gene.__dict__)


        if not gene.transcripts or not gene.cdss or not gene.coordinates:
            # no transcript was found for gene
            continue

        seq = ''
        # assess sequence with considering the phase attribute and without it
        seqs = {'phased': '','non-phased': ''}
        for coordinate in gene.coordinates:
            # on reverse strand CDS needs to be reverse complemented before concatenated
            if gene.strand:
                seqs['non-phased'] += str(Seq.Seq(contig_seq[coordinate[0]:coordinate[1]]).reverse_complement())
            else:
                seqs['non-phased'] += str(Seq.Seq(contig_seq[coordinate[0]:coordinate[1]]))
        for coordinate in gene.phased_coordinates:
            # on reverse strand CDS needs to be reverse complemented before concatenated
            if gene.strand:
                seqs['phased'] += str(Seq.Seq(contig_seq[coordinate[0]:coordinate[1]]).reverse_complement())
            else:
                seqs['phased'] += str(Seq.Seq(contig_seq[coordinate[0]:coordinate[1]]))

        if len(seqs['phased'])%3 != 0:
            outta_frame['p'] += 1
            seqs['phased'] = seqs['phased']+("N"*(3-(len(seqs['phased'])%3)))
        if len(seqs['non-phased'])%3 != 0:
            outta_frame['np'] += 1
            seqs['non-phased'] = seqs['non-phased']+("N"*(3-(len(seqs['non-phased'])%3)))

        protein = str(Seq.Seq(seqs['non-phased']).translate(table=int(gene.transl_table)))
        phased_protein = str(Seq.Seq(seqs['phased']).translate(table=int(gene.transl_table)))

        gene.prot_seqs['phased'] = phased_protein
        gene.prot_seqs['non-phased'] = protein



        if "*" in phased_protein[-1]:
            stop_codons['p'] += 1
        if "*" in protein[:-1]:
            stop_codons['np'] += 1


    return stop_codons, outta_frame


def decide_phasing(cfg, stop_codons, outta_frame):
    """Identify if information of phase field should be used"""
    if cfg.use_phase_info == 'auto':
        logging.debug(f"Proteins with phase information: "
                      f"\n\tinternal stop codons: {stop_codons.get('p')}"
                      f"\n\tproteins out of frame: {outta_frame.get('p')}")
        logging.debug(f"Proteins without phase information: "
                      f"\n\tinternal stop codons: {stop_codons.get('np')}"
                      f"\n\tproteins out of frame: {outta_frame.get('np')}")
        phasing = stop_codons.get('p') + outta_frame.get('p')
        non_phasing = stop_codons.get('np') + outta_frame.get('np')

        if phasing > non_phasing:
            return False
        if phasing < non_phasing:
            return True
        else: # both are equal
            # weigh internal stop codons stronger
            if stop_codons.get('p') > stop_codons.get('np'):
                return False
            if stop_codons.get('p') < stop_codons.get('np'):
                return True
            else:
                # per default use phase information if both score equally good
                return True

    if cfg.use_phase_info:
        return True
    if not cfg.use_phase_info:
        return False

def write_seqs(cfg, contigs, features, use_phase_info):
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

    with open(cfg.proteins_path, "w") as proteins_file:

        for gene_lists in contigs.values():
            for gene_id in gene_lists:
                gene = features.get(gene_id)
                # print(gene.__dict__)

                proteins_file.write(">"+gene_id+"\n")

                if use_phase_info:
                    protein = gene.prot_seqs.get('phased')
                else:
                    protein = gene.prot_seqs.get('non-phased')

                proteins_file.write(protein)
                proteins_file.write("\n")


def extract_seq(cfg, contigs, features):
    """Retrieve sequence for contig from FASTA and write protein FASTA.

    Read assembly FASTA, retrieve scaffold sequence and call set_seqs()
    to write AA sequences for genes to file

    Args:
      cfg(obj): Config class object of config parameters
      contigs(dict): {scaffold ID : [gene IDs]}
      features(dict): {GFF feature ID : Feature object}
    """

    stop_codons = {'p': 0, 'np': 0}
    outta_frame = {'p': 0, 'np': 0}

    current_contig = ""

    with open(cfg.fasta_path, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                # False if no genes in contig dict (geneless contigs)
                if contigs.get(current_contig):
                    stop_codons, outta_frame = set_seqs(contigs, features, current_contig, contig_seq,stop_codons, outta_frame)
                # retrieve ID from current contig
                current_contig = line.strip().lstrip(">").split()[0]
                contig_seq = ''
            else:
                contig_seq += line.strip()
        # add genes from last contig
        if contigs.get(current_contig):
            stop_codons, outta_frame =  set_seqs(contigs, features, current_contig, contig_seq,stop_codons, outta_frame)

    return stop_codons, outta_frame

def log_seqs_info(use_phase_info, stop_codons, outta_frame):

    if use_phase_info:
        internal_stops = stop_codons.get('p')
        outta_frames = outta_frame.get('p')
    else:
        internal_stops = stop_codons.get('np')
        outta_frames = outta_frame.get('np')

    errors = False
    if internal_stops != 0:
        logging.warning(f'{internal_stops} proteins with '
                        f'internal stop codon(s) identified')
        errors = True
    if outta_frames != 0:
        logging.warning(f'{outta_frames} proteins with partial codon(s) '
                        f'identified (length of sequence no multiple of three).'
                        f' Trailing Ns were added')
        errors = True
    if errors:
        logging.warning(
            'This may be a problem with your GFF. If you have a precomputed '
            'protein FASTA file available, using this may yield a better '
            'taxonomic assignment (especially if the number(s) is/are high).')


def generate_fasta(cfg):
    """Generate FASTA file of protein sequences.

    Extract coordinates for genes and transcripts and their relation,
    determine the longest transcript for each gene,
    extract the nucleotide sequence and translate to amino acid sequence

    Args:
      cfg(obj): Config class object of config parameters
    """

    contigs, features = get_cds_coordinates(cfg)
    get_longest_transcript(contigs, features)
    stop_codons, outta_frame = extract_seq(cfg, contigs, features)
    use_phase_info = decide_phasing(cfg, stop_codons, outta_frame)
    write_seqs(cfg, contigs, features, use_phase_info)
    log_seqs_info(use_phase_info, stop_codons, outta_frame)



def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkData.cfg2obj(config_path)

    generate_fasta(cfg)
    process_coverage(cfg)


if __name__ == '__main__':
    main()
