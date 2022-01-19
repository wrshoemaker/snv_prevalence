#!/bin/bash
#$ -N predict_prevalence
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/prevalence_error_test
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/prevalence_output_test
#$ -l h_data=24G
#$ -l time=120:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


#species=$1

python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --make_alleles_dict all

#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --make_alleles_dict Alistipes_finegoldii_56071

#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --species Alistipes_finegoldii_56071

#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --species Bacteroides_xylanisolvens_57185


#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --make_frequency_dict_mapgd_non_zero


#qsub -N Bacteroides_xylanisolvens_57185 /u/home/w/wrshoema/project-ngarud/snv_prevalence/scripts/qsub_prevalence/qsub_predicted_prevalence.sh Bacteroides_xylanisolvens_57185


# qrsh -l h_rt=3:00:00,h_data=24G
