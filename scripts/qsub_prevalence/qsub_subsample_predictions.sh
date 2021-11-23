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
module load python/2.7.15


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence.py --estimate_pi

#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence.py --subsample_predictions


# qrsh -l h_rt=3:00:00,h_data=24G
