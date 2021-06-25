#!/bin/bash
#$ -N plot_for_meeting
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_for_meeting_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_for_meeting_output
#$ -l h_data=128G
#$ -l time=24:00:00
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/3.6.1

python3 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/plot_for_meeting.py
