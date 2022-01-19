#!/bin/bash
#$ -N simulation
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/simulation_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/simulation_output
#$ -l h_data=6G
#$ -l time=16:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


f_mean=$1
beta=$2

python /u/project/ngarud/wrshoema/snv_prevalence/scripts/simulate_gamma_estimator.py --f_mean ${f_mean} --beta ${beta}



#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/simulate_gamma_estimator.py  --f_mean 0.1 --beta 1 --n_iter 1


# qrsh -l h_rt=3:00:00,h_data=24G
