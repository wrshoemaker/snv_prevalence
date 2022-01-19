#!/bin/bash
#$ -N Barnesiella_intestinihominis_62208
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/Barnesiella_intestinihominis_62208_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/Barnesiella_intestinihominis_62208_output
#$ -l h_data=24G
#$ -l time=120:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --species Barnesiella_intestinihominis_62208


# qrsh -l h_rt=3:00:00,h_data=24G
