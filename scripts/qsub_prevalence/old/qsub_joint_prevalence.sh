#!/bin/bash
#$ -N joint_prevalence
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_joint_prevalence_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_joint_prevalence_output
#$ -l h_data=36G
#$ -l time=72:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13


python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_joint_prevalence.py
