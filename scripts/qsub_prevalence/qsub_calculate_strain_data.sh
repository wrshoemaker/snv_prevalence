#!/bin/bash
#$ -N calculate_ld
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_prevalence_summary
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_prevalence_summary
#$ -l h_data=8G
#$ -l time=4:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
#module load python/3.6.1
module load python/2.7.13


python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_prevalence_summary.py
