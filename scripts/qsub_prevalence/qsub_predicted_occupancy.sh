#!/bin/bash
#$ -N predict_prevalence
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_occupancy_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_occupancy_output
#$ -l h_data=24G
#$ -l time=72:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_occupancy.py
