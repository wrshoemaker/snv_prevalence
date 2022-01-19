#!/bin/bash
#$ -N Bacteroides_caccae_53434
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/Bacteroides_caccae_53434_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/Bacteroides_caccae_53434_output
#$ -l h_data=24G
#$ -l time=120:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --species Bacteroides_caccae_53434


# qrsh -l h_rt=3:00:00,h_data=24G
