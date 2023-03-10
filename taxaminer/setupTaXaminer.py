#!/usr/bin/env python


"""setup external dependencies for taXaminer

"""

__author__ = "Freya Arthen"
__version__ = "0.6.0"

import argparse
import os
import pathlib
import stat
import tarfile
import zipfile
import requests
import shutil
import subprocess
import sys
import glob


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
    print("checking diamond")
    check_diamond = subprocess.run([f"{cmd_dict.get('diamond')} --version"],
                                     shell=True, capture_output=True)
    if check_diamond.returncode != 0:
        print("Installation of diamond failed. Please try again "
                        "by running\n\nconda install -c bioconda diamond\n")

    # check Krona tools
    print("checking krona tools")
    check_krona = subprocess.run([f"{cmd_dict.get('krona')}"],
                                   shell=True, capture_output=True)
    if check_krona.returncode != 0:
        print("Installation of Krona failed. Please try again "
                        "by running\n\nconda install -c bioconda krona\n")


def setup_db(outPath, cmd_dict):
    pathlib.Path(outPath).mkdir(parents=True, exist_ok=True)
    print("downloading nr")
    download_file(f"{outPath}/nr.gz", "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz")
    print("downloading taxdmp")
    download_file(f"{outPath}/taxdump.tar.gz",
                  "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
    print("downloading accession")
    download_file(f"{outPath}/prot.accession2taxid.FULL.gz",
                  "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz")

    print("unpack")
    open_tarfile(f"{outPath}/taxdump.tar.gz", outPath, "gz", ["names.dmp", "nodes.dmp", "merged.dmp"])

    print("preparing database")
    index_db = (f'{cmd_dict.get("diamond")} makedb --in "{outPath}/nr.gz" -d "{outPath}/nr" \
        --taxonmap "{outPath}/prot.accession2taxid.FULL.gz" \
        --taxonnodes "{outPath}/nodes.dmp" \
        --taxonnames "{outPath}/names.dmp" \
        && rm -rf nr.gz prot.accession2taxid.gz taxdump.tar.gz')
    out_index_db = subprocess.run([index_db], shell=True, capture_output=True)
    if out_index_db.returncode != 0:
        print(
            f"creation of diamond database failed:\n{out_index_db}")
        print("Error message:\n" + out_index_db.stderr.decode())
        sys.exit()

def setup_krona_taxonomy(path):
    print("setting up taxonomy db for krona tools")
    setup_krona = subprocess.run([f"{path}"],
                                 shell=True, capture_output=True)
    if setup_krona.returncode != 0:
        print("Setup of Krona database failed. Error message:"
              f"{setup_krona.stderr.decode()}")


def setup_conda():
    # check conda
    check_conda = subprocess.run(["conda --version"], shell=True, capture_output=True)
    if check_conda.returncode != 0:
        print("Could not find anaconda installation.")
        sys.exit()

    cmd_dict = {}

    # install diamond
    print("installing diamond")
    install_diamond = subprocess.run(["conda install -c bioconda diamond -y"],
                                     shell=True, capture_output=True)
    if install_diamond.returncode != 0:
        print("Installation of diamond failed. Please try again "
                        "by running\n\nconda install -c bioconda diamond\n")
    cmd_dict["diamond"] = "diamond"

    # install Krona tools
    print("installing krona tools")
    install_krona = subprocess.run(["conda install -c bioconda krona -y"],
                                   shell=True, capture_output=True)
    if install_krona.returncode != 0:
        print("Installation of Krona failed. Please try again "
                        "by running\n\nconda install -c bioconda krona\n")
    cmd_dict["krona"] = "ktImportTaxonomy"
    ## initializing taxonomy for Krona
    krona_run_path = subprocess.run(["which ktImportTaxonomy"], shell=True,
                                capture_output=True).stdout.decode()
    krona_update_path = krona_run_path.replace('bin/ktImportTaxonomy',
                                               'opt/krona/updateTaxonomy.sh')
    setup_krona_taxonomy(krona_update_path)


    return cmd_dict


def setup_locally(tool_path):
    # make dir for tools
    toolPath = pathlib.Path(tool_path).absolute()
    pathlib.Path(toolPath).mkdir(parents=True, exist_ok=True)
    os.chdir(toolPath)
    toolPath = str(toolPath) + "/"

    print(toolPath)

    cmd_dict = {}

    # install diamond
    print("installing diamond")
    dmd_version = "2.0.15"
    download_file(f"diamond-linux64.tar.gz",
                  "https://github.com/bbuchfink/diamond/releases/download/"
                  f"v{dmd_version}/diamond-linux64.tar.gz")
    open_tarfile(f"diamond-linux64.tar.gz", toolPath, "gz")
    cmd_dict["diamond"] = f"{toolPath}diamond"

    # install Krona tools
    print("installing krona tools")
    k_version = "2.8.1"
    download_file(f"KronaTools-{k_version}.tar",
                  "https://github.com/marbl/Krona/releases/download/"
                  f"v{k_version}/KronaTools-{k_version}.tar")
    open_tarfile(f"KronaTools-{k_version}.tar", toolPath, "")
    os.chdir(pathlib.Path(f"KronaTools-{k_version}/"))
    install_krona = subprocess.run([f"./install.pl --prefix .;"],
                                   shell=True, capture_output=True)
    os.chdir(toolPath)
    if install_krona.returncode != 0:
        print("Installation of Krona failed:\n"
                        f"{install_krona.stderr.decode()}")
    else:
        os.symlink(f"KronaTools-{k_version}/bin/ktImportTaxonomy",
                   "ktImportTaxonomy")
    cmd_dict["krona"] = f"{toolPath}ktImportTaxonomy"
    ## initializing taxonomy for Krona
    krona_update_path = f"{toolPath}KronaTools-{k_version}/updateTaxonomy.sh"
    setup_krona_taxonomy(krona_update_path)

    # remove downloads
    os.remove("diamond-linux64.tar.gz")
    os.remove(f"KronaTools-{k_version}.tar")

    return cmd_dict

def main():
    """  """
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
                          action="store_true", default=False)
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
            sys.exit("No config file found. Please run taxaminer.setup with option '-db'")
        else:
            print(taxaminer_paths_save.get("data_path"))
            sys.exit()
    if args.getSourcepath:
        print(taxaminer_path)
        sys.exit()
    if args.getToolpath:
        if not taxaminer_paths_save:
            sys.exit("No config file found. Please run taxaminer.setup")
        else:
            print(taxaminer_paths_save.get("tool_path"))
            sys.exit()

    # sanity check of user input
    if (not args.conda and not args.toolPath) and (args.database and not args.dataPath):
        sys.exit("Path for the installation not given!\nPlease use option "
                "'-o' to specify the path where tools will be installed\n"
                "and '-d' for path where database will be stored")

    if args.conda:
        print("setting up dependencies via conda")
        cmd_dict = setup_conda()
        tool_path = "conda installation"
        check_executable(cmd_dict)
    elif args.toolPath:
        print("setting up dependencies locally")
        if not args.toolPath and not taxaminer_paths_save:
            sys.exit("No path for the installation given. Please use option"
                     "'-o' to specify the path where tools will be installed")
        if args.toolPath:
            tool_path = os.path.abspath(args.toolPath)
        else:
            tool_path = taxaminer_paths_save.get("tool_path")
            # TODO: delete existing versions of databases?
            # shutil.rmtree(tool_path)
        cmd_dict = setup_locally(tool_path)
        check_executable(cmd_dict)

    #TODO: Krona: add update_taxonomy

    if args.database:
        print("setting up database")
        if not tool_path:
            sys.exit("Tools for database creation not installed. Please run "
                     "taxaminer.setup with either option '-o' or '--conda'.")
        if not args.dataPath and not taxaminer_paths_save:
            sys.exit("No path for database given. Please use option"
                     "'-d' to specify the database will be stored.")
        if args.dataPath:
            data_path = os.path.abspath(args.dataPath)
        else:
            data_path = taxaminer_paths_save.get("data_path")
        # TODO: use tools either via conda or via local installation
        setup_db(args.dataPath, cmd_dict)
    elif args.dataPath:
        # database already precomputed -> adjust the path
        data_path = os.path.abspath(args.dataPath)

    taxaminer_paths = {"data_path": data_path,
                        "tool_path": tool_path,
                       "cmds": cmd_dict}

    with open(taxaminer_path+"/pathconfig.txt", "w") as config:
            config.write(str(taxaminer_paths))
            config.close()


if __name__ == "__main__":
    main()
