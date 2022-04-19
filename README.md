[![DOI:10.1101/2021.10.05.463161](http://img.shields.io/badge/DOI-10.1101/2021.10.05.463161-B31B1B.svg)](https://doi.org/10.1101/2022.04.07.487434)


# snv_prevalence

Repository for code associated with the preprint:

[A macroecological perspective on genetic diversity in the human gut microbiome](https://doi.org/10.1101/2022.04.07.487434).



## Dependencies
All code was written in Python 2.7. Details of the conda environment can be found in `environment.yml`. Set up the repo under a folder named Github: `~/GitHub/snv_prevalence/`.


## Running the analyses

All the commands to reprocess the data are called in `run_everything.sh`. This script contains commands for all data processing and analyses, including MIDAS. You will need to rework the paths in the scripts listed in `run_everything.sh` to get it working on your machine, including moving all HMP fastq files into `~/GitHub/snv_prevalence/data/SRA_files/`. All bash scripts are written to be run on a cluster. I do not recommend running them on your local machine.


Python scripts in `run_everything.sh` that are marked for figure generation can be run on a local machine. These scripts require the processed data.
