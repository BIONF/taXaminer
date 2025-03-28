# taXaminer

taXaminer - examine the taxonomic diversity in genome assemblies. Designed to detect and differentiate contamination and horizontal gene transfer. 

taXaminer combines a reference-free and an alignment-based approach to detect and differentiate contamination and horizontal gene transfer in genome assemblies. It uses a total of 16 intrinsic features to describe the gene set. Among these are the read coverage, sequence composition, gene length and the size of the scaffold it is annotated on (see details [here](https://github.com/BIONF/taXaminer/wiki/Additional-information#pca-variables)). To identify genes which discern from the average, a Principal Component Analysis is used to cluster genes with similar features. The taxonomic assignment targets at identifying the true taxon of origin for each gene. It is based on their protein sequence to reduce the need of having the exact reference in the database. 

The results can be interactively explored in the accompanying [dashboard](https://github.com/BIONF/taxaminer-dashboard).

# Table of Contents
* [Installation](#installation)
* [Usage](#usage)
* [Bugs](#bugs)
* [Contributors](#contributors)
* [License](#license)
* [Contact](#contact)

# Installation

To install taXaminer, use the python package installer pip. 
Note: taXaminer is as of yet not published at pypi, thus you need to download this repository and provide pip with the link to the directory for installation.
```
git clone https://github.com/BIONF/taXaminer.git
pip install ./taXaminer
```

To install the additional dependencies, use the setup function included in taXaminer. You can install the tools either via conda or locally in a specified directory.

Using conda (installs into the currently active environment):
```
taxaminer.setup --conda
```
In a local directory:
```
taxaminer.setup -o </path/to/tool/directory/>
```
To download and build the database, use:
```
taxaminer.setup -db -d </path/to/database/directory/>
```
Use the following command to use an existing database.
```
taxaminer.setup -d </path/to/existing_database/directory/>
```

# Usage
1. Create a configuration file using the following template and adapt it to fit your data.
```
fasta_path: "path/to/assembly.fasta" # path to assembly FASTA
gff_path: "path/to/assembly.gff" # path to annotation in GFF3 format
output_path: "path/to/output_directory/" # directory to save results in
taxon_id: "<NCBI taxon ID>" # NCBI Taxon ID of query species
```
2. To include coverage information, add the path to a sorted bam file (this is optional). Otherwise, omit this parameter from the configuration file.
```
bam_path_1: "path/to/mapping.bam" # path to BAM file
```
* Note: When using multiple coverage sets, duplicate the parameter and increase the number in the suffix


To run taXaminer, call it with the path to the config file, like so:
```
taxaminer.run <config.yml>
```

For details on additional options see [Configuration parameters](https://github.com/BIONF/taXaminer/wiki/Configuration-parameters). 

# Bugs
Any bug reports, comments or suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/taXaminer/issues/new) or reach out via [email](mailto:f.arthen@bio.uni-frankfurt.de).

# Contributors
* [Freya Arthen](https://github.com/fdarthen)
* Simonida Zehr

# License
*taXaminer* is released under [MIT license](https://github.com/BIONF/taXaminer/blob/master/LICENSE).

# Contact
Please contact us via [email](mailto:f.arthen@bio.uni-frankfurt.de).
