#!/bin/bash
#$ -N get_coverage_alt_allele_mapgd
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/scripts/get_coverage_alt_allele_mapgd_error
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/scripts/get_coverage_alt_allele_mapgd_output
#$ -l h_data=24G
#$ -l time=48:00:00
#$ -l highp
#$ -m bea


. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15


python /u/project/ngarud/wrshoema/snv_prevalence/scripts/get_coverage_alt_allele_mapgd.py
