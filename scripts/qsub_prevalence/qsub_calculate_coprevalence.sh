#!/bin/bash
#$ -N calculate_coprevalence
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence_output
#$ -l h_data=8G
#$ -l time=240:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13


#python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence.py

python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence_f0.py
