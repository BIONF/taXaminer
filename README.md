# taXaminer

taXaminer enables the interactive exploration of taxonomic footprints in gene sets. The specific goal is to detect and differentiate contamination and horizontal gene transfer.

Besides the taxonomic assignment of genes, taXaminer uses a total of 16 further indicators to faciliate this. Among these indicators are read coverage, sequence composition, gene length and position of genes within their scaffold. To identify genes which deviate from the mean set of genes, a principal component analysis (PCA) is used as it condenses data to fewer dimensions. Genes with similar values for certain variables are thereby clustered together, so that deviations are made visible. The results can be interactively examined in a 3D scatterplot, where the dot position respresents a combination of coverage, sequence composition and spatial information provided by the PCA and the color the taxonomic assignment.

# Table of Contents
* [Installation](#installation)
* [Usage](#usage)
* [Bugs](#bugs)
* [Contributors](#contributors)
* [License](#license)
* [Contact](#contact)

# Installation

To install taXaminer, use the python package installer pip.
```
pip install taxaminer
```

To install additional dependencies, use the setup function included in taXaminer. You can install the tools either via conda or locally in a specified directory.

Using conda:
```
taxaminer.setup --conda
```
In a local directory:
```
taxaminer.setup --o </path/to/tool/directory/>
```
To download and build the database, use:
```
taxaminer.setup --db -d </path/to/database/directory/>
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
2. To include coverage information, add one of the following parameters according to your data (this is optional).
```
pbc_path_1: "path/to/pbc.txt" # path to PBC file; omit to use default location in output directory
bam_path_1: "path/to/mapping.bam" # path to BAM file; omit to use default location in output directory
read_paths_1: ["path/to/read_file_1.fa","path/to/read_file_2.fa"] # path to read file(s)
```
* Note: When using multiple coverage sets, duplicate the parameter you need and increase the number in the suffix


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
