#!/bin/bash
#$ -N plot_for_meeting
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_for_meeting_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_for_meeting_output
#$ -l h_data=24G
#$ -l time=48:00:00
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13


python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_dissimilarity_overlap.py
