#!/bin/bash
#$ -N get_maf_all_species
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/get_maf_all_species_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/get_maf_all_species_output
#$ -l h_data=8G
#$ -l time=8:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13

python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/get_maf_all_species.py
