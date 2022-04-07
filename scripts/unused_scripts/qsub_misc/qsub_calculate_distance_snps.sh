#!/bin/bash
#$ -N distance_snps
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_distance_snps_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_distance_snps_output
#$ -l h_data=28G
#$ -l time=48:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13


python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_distance_snps.py
