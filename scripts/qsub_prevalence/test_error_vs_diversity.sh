#!/bin/bash
#$ -N predict_prevalence
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/error_vs_diversity_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/error_vs_diversity_output
#$ -l h_data=24G
#$ -l time=72:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/plot_error_vs_diversity.py




#qrsh -l h_rt=1:00:00,h_data=16G
