#!/bin/bash
#$ -N Eubacterium_rectale_56927
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/Eubacterium_rectale_56927_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/Eubacterium_rectale_56927_output
#$ -l h_data=24G
#$ -l time=120:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --species Eubacterium_rectale_56927


# qrsh -l h_rt=3:00:00,h_data=24G
