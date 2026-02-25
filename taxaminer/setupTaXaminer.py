#!/usr/bin/env python


"""setup external dependencies for taXaminer

"""

__author__ = "Freya Arthen"

import argparse
import gzip
import os
import pathlib
import tarfile
import requests
import subprocess
import sys
import logging


def download_file(local_filename, url):
    #local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, "wb") as file:
            for chunk in r.iter_content(chunk_size=8192):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                file.write(chunk)
    return local_filename


def open_tarfile(file, destination, compression, files=[]):
    tar = tarfile.open(file, f"r:{compression}")
    if files:
        for file in files:
            tar.extract(file, destination)
    else:
        tar.extractall(destination)
    tar.close()


def write_tool_cmds(taxaminer_path, cmd_dict):

    with open(f"{taxaminer_path}/pathconfig.txt", "a") as config_file:
        config_file.write(str(cmd_dict))


def check_executable(cmd_dict):

    #TODO: check what this should print out

    # check diamond
    logging.info("checking diamond")
    check_diamond = subprocess.run([f"{cmd_dict.get('diamond')} --version"],
                                     shell=True, capture_output=True)
    if check_diamond.returncode != 0:
        logging.error("Installation of diamond failed. Please try again "
                        "by running\n\nconda install -c bioconda diamond\n")


def prepare_nr(outPath, db_name):
    logging.info(f">> downloading {db_name}")
    download_file(f"{outPath}/{db_name}.gz",
                  f"https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/{db_name}.gz")

    logging.info(">> downloading accession mapping file")
    download_file(f"{outPath}/prot2taxid.tsv.gz",
                  "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz")

    return f"{outPath}/{db_name}.gz"


def prepare_uniref(outPath, db_name):
    logging.info(f">> downloading {db_name}")
    download_file(f"{outPath}/{db_name}.gz",
                  f"https://ftp.uniprot.org/pub/databases/uniprot/uniref/{db_name}/{db_name}.fasta.gz")
    logging.info(">> downloading accession mapping file")
    download_file(f"{outPath}/idmapping_selected.tab.gz",
                  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz")

    # Columns idmapping_selected.tab
    # # 0.     UniProtKB - AC
    # # 1.     UniProtKB - ID
    # # 2.     GeneID(EntrezGene)
    # # 3.     RefSeq
    # # 4.     GI
    # # 5.     PDB
    # # 6.     GO
    # # 7.     UniRef100
    # # 8.     UniRef90
    # # 9.     UniRef50
    # # 10.    UniParc
    # # 11.    PIR
    # # 12.    NCBI - taxon
    # # 13.    MIM
    # # 14.    UniGene
    # # 15.    PubMed
    # # 16.    EMBL
    # # 17.    EMBL - CDS
    # # 18.    Ensembl
    # # 19.    Ensembl_TRS
    # # 20.    Ensembl_PRO
    # # 21.    Additional PubMed

    # read taxids that are in nodes files
    nodes_taxids = set()
    with open(f"{outPath}/nodes.dmp", 'r') as nodes_file:
        for node in nodes_file:
            nodes_taxids.add(node.split()[0])
    old2new = {}
    with open(f"{outPath}/merged.dmp", 'r') as nodes_file:
        for match in nodes_file:
            old2new[match.split()[0]] = match.split()[1]

    # select correct column for respective uniref version
    if db_name == 'Uniref100':
        col_num = 7
    elif db_name == 'Uniref90':
        col_num = 8
    else:
        col_num = 9

    mapping_dict = {}
    logging.info(">> creating protein to taxon mapping file")
    with gzip.open(f"{outPath}/idmapping_selected.tab.gz", 'rb') as in_mapping, \
        gzip.open(f"{outPath}/prot2taxid.tsv.gz", 'wb', compresslevel=6) as out_mapping:
        out_mapping.write(f"accession.version\ttaxid\n".encode('utf-8'))
        for line in in_mapping:
            spline = line.decode('utf-8').strip().split('\t')
            if spline[col_num] and spline[12]:
                if spline[12] in nodes_taxids:
                    out_mapping.write(f"{spline[0]}\t{spline[12]}\n".encode('utf-8'))
                    if spline[col_num] in mapping_dict.keys():
                        mapping_dict[spline[col_num]].append(spline[0])
                    else:
                        mapping_dict[spline[col_num]] = [spline[0]]
                elif spline[col_num] in old2new.keys():
                    new_taxid = old2new.get(spline[col_num])
                    out_mapping.write(
                        f"{spline[0]}\t{new_taxid}\n".encode('utf-8'))
                    if spline[col_num] in mapping_dict.keys():
                        mapping_dict[spline[col_num]].append(spline[0])
                    else:
                        mapping_dict[spline[col_num]] = [spline[0]]

    logging.info(">> renaming fasta headers in database file")
    with gzip.open(f"{outPath}/{db_name}.gz", 'rb') as in_db_file, \
        gzip.open(f"{outPath}/db.gz", 'wb', compresslevel=6) as out_db_file:
        for line in in_db_file:
            if line.decode('utf-8').startswith('>'):
                uniref_id = line.decode('utf-8').split()[0][1:]
                if uniref_id in mapping_dict.keys():
                    joined_accessions = '\x01'.join(mapping_dict.get(uniref_id))
                    out_db_file.write(f">{joined_accessions}\n".encode('utf-8'))
                else:
                    prefix_stripped_line = '_'.join(line.decode('utf-8').split('_')[1:])
                    out_db_file.write(f">{prefix_stripped_line}".encode('utf-8'))
            else:
                out_db_file.write(line)

    os.remove(f"{outPath}/idmapping_selected.tab.gz")

    return f"{outPath}/db.gz"

def setup_db(db_name, outPath, cmd_dict):
    pathlib.Path(outPath).mkdir(parents=True, exist_ok=True)
    logging.info("> downloading files")
    logging.info(">> downloading taxdmp")
    download_file(f"{outPath}/taxdump.tar.gz",
                  "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
    logging.info("> unpacking taxdmp")
    open_tarfile(f"{outPath}/taxdump.tar.gz", outPath, "gz",
                 ["names.dmp", "nodes.dmp", "merged.dmp"])

    if db_name == 'nr' or db_name == 'swissprot':
        db_path = prepare_nr(outPath, db_name)
    elif "uniref" in db_name:
        db_path = prepare_uniref(outPath, db_name)


    logging.info("> preparing database")
    index_db = (f'{cmd_dict.get("diamond")} makedb -v --in "{db_path}" '
                f'-d "{outPath}/db" '
                f'--taxonmap "{outPath}/prot2taxid.tsv.gz" '
                f'--taxonnodes "{outPath}/nodes.dmp" '
                f'--taxonnames "{outPath}/names.dmp" ')
                #f'&& rm -rf {outPath}/{db_name}.gz {outPath}/db.gz {outPath}/prot2taxid.tsv.gz {outPath}/taxdump.tar.gz')
    out_index_db = subprocess.run([index_db], shell=True, capture_output=True)
    if out_index_db.returncode != 0:
        logging.info(f"creation of diamond database failed:\n{out_index_db}")
        logging.info("Error message:\n" + out_index_db.stderr.decode())
        sys.exit()
    else:
        logging.info(f">> creation of diamond database successful!")
        logging.info(
            f"Setup log:\n{out_index_db}")
        logging.info("Stdout:\n" + out_index_db.stderr.decode())


def setup_conda():
    # check conda
    check_conda = subprocess.run(["conda --version"], shell=True, capture_output=True)
    if check_conda.returncode != 0:
        logging.error("Could not find anaconda installation.")
        sys.exit()

    cmd_dict = {}

    # install diamond
    logging.info("installing diamond")
    install_diamond = subprocess.run(['conda install -c bioconda -c conda-forge "diamond>=2.1.7" -y'],
                                     shell=True, capture_output=True)
    if install_diamond.returncode != 0:
        logging.error("Installation of diamond failed. Please try again "
                        'by running\n\nconda install -c bioconda "diamond>=2.1.7"\n')
    cmd_dict["diamond"] = "diamond"

    logging.info("installing chrome for kaleido")    
    import kaleido
    kaleido.get_chrome_sync()

    return cmd_dict


def setup_locally(tool_path):
    # make dir for tools
    toolPath = pathlib.Path(tool_path).absolute()
    pathlib.Path(toolPath).mkdir(parents=True, exist_ok=True)
    os.chdir(toolPath)
    toolPath = str(toolPath) + "/"

    cmd_dict = {}

    # install diamond
    logging.info("installing diamond")
    dmd_version = "2.0.15"
    download_file(f"diamond-linux64.tar.gz",
                  "https://github.com/bbuchfink/diamond/releases/download/"
                  f"v{dmd_version}/diamond-linux64.tar.gz")
    open_tarfile(f"diamond-linux64.tar.gz", toolPath, "gz")
    cmd_dict["diamond"] = f"{toolPath}diamond"

    # remove downloads
    os.remove("diamond-linux64.tar.gz")

    logging.info("installing chrome for kaleido")    
    import kaleido
    kaleido.get_chrome_sync()

    return cmd_dict

def main():
    """  """

    logging.basicConfig(level=logging.INFO,  # Set the logging level to INFO
                        format='%(levelname)s: %(message)s')

    # Create a logger object
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    optional.add_argument("-o", "--toolPath",
                          help="Path for installation of external tools",
                          action="store", default="")
    optional.add_argument("--conda", help="Installing external tools within a conda env",
                          action="store_true", default=False)
    optional.add_argument("-db", "--database",
                          help="Download/Update database for taXaminer",
                          choices=['nr', 'swissprot', 'uniref90', 'uniref100', 'uniref50'])
    optional.add_argument("-d", "--dataPath",
                          help="Path to taXaminer database",
                          action="store", default="")
    optional.add_argument("--getSourcepath",
                          help="Get path to installion of taXaminer",
                          action="store_true", default=False)
    optional.add_argument("--getDatapath",
                          help="Get taXaminer path to database",
                          action="store_true", default=False)
    optional.add_argument("--getToolpath",
                          help="Get taXaminer path to installed tools",
                          action="store_true", default=False)

    args = parser.parse_args()

    taxaminer_path = os.path.realpath(__file__).replace("setupTaXaminer.py","")
    # check if setting have already been made
    taxaminer_config = taxaminer_path+"/pathconfig.txt"
    if os.path.exists(taxaminer_config):
        with open(taxaminer_config) as f:
            taxaminer_paths_save = eval(f.readline().strip())
            data_path = taxaminer_paths_save.get("data_path")
            tool_path = taxaminer_paths_save.get("tool_path")
            if taxaminer_paths_save.get("cmds"):
                cmd_dict = dict(taxaminer_paths_save.get("cmds"))
            else:
                cmd_dict = None
    else:
        taxaminer_paths_save = None
        data_path = None
        tool_path = None
        cmd_dict = None

    # print requested paths
    if args.getDatapath:
        if not taxaminer_paths_save:
            logger.error("No config file found. Please run taxaminer.setup with option '-db'")
            sys.exit()
        else:
            logging.info(taxaminer_paths_save.get("data_path"))
            sys.exit()
    if args.getSourcepath:
        logging.info(taxaminer_path)
        sys.exit()
    if args.getToolpath:
        if not taxaminer_paths_save:
            logger.error(
                "No config file found. Please run taxaminer.setup")
            sys.exit()
        else:
            logging.info(taxaminer_paths_save.get("tool_path"))
            sys.exit()

    # sanity check of user input
    if (not args.conda and not args.toolPath) and (args.database and not args.dataPath):
        logger.error("Path for the installation not given!\nPlease use option "
                    "'-o' to specify the path where tools will be installed\n"
                    "and '-d' for path where database will be stored")
        sys.exit()

    if args.conda:
        logging.info("setting up dependencies via conda")
        cmd_dict = setup_conda()
        tool_path = "conda installation"
        check_executable(cmd_dict)
    elif args.toolPath:
        logging.info("setting up dependencies locally")
        if not args.toolPath and not taxaminer_paths_save:
            logger.error("No path for the installation given. Please use option"
                     "'-o' to specify the path where tools will be installed")
            sys.exit()
        if args.toolPath:
            tool_path = os.path.abspath(args.toolPath)
        else:
            tool_path = taxaminer_paths_save.get("tool_path")
            # TODO: delete existing versions of databases?
            # shutil.rmtree(tool_path)
        cmd_dict = setup_locally(tool_path)
        check_executable(cmd_dict)

    if args.database:
        logging.info("setting up database")
        if not tool_path:
            logger.error("Tools for database creation not installed. Please run "
                     "taxaminer.setup with either option '-o' or '--conda'.")
            sys.exit()
        if not args.dataPath and not taxaminer_paths_save:
            logger.error("No path for database given. Please use option"
                     "'-d' to specify the database will be stored.")
            sys.exit()
        if args.dataPath:
            data_path = os.path.abspath(args.dataPath)
        else:
            data_path = taxaminer_paths_save.get("data_path")
        # TODO: use tools either via conda or via local installation
        setup_db(args.database, args.dataPath, cmd_dict)
    elif args.dataPath:
        # database already precomputed -> adjust the path
        logger.info(f"Setting database path to \n\t{os.path.abspath(args.dataPath)}")
        data_path = os.path.abspath(args.dataPath)

    taxaminer_paths = {"data_path": data_path,
                        "tool_path": tool_path,
                       "cmds": cmd_dict}

    with open(taxaminer_path+"/pathconfig.txt", "w") as config:
            config.write(str(taxaminer_paths))
            config.close()

    logging.info("Done!")


if __name__ == "__main__":
    main()
