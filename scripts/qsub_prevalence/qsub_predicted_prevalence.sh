#!/bin/bash
#$ -N predict_prevalence
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/prevalence_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/prevalence_output
#$ -l h_data=24G
#$ -l time=72:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


species=$1

python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --species ${species}

#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --make_alleles_dict ${species}

#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --make_alleles_dict ${species}



#python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --species Alistipes_finegoldii_56071


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/calculate_predicted_prevalence_mapgd.py --make_frequency_dict_mapgd_non_zero

# qrsh -l h_rt=3:00:00,h_data=24G
