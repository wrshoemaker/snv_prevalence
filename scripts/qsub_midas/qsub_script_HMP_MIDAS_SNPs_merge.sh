#!/bin/bash
#$ -N snps_merge
#$ -e /u/project/ngarud/Garud_lab/snv_prevalence/scripts/snps_error
#$ -o /u/project/ngarud/Garud_lab/snv_prevalence/scripts/snps_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=120:00:00
#$ -l h_data=80G
#$ -l highp


. /u/local/Modules/default/init/modules.sh


module unload python
module load anaconda/python2-4.2

source activate midas


export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2


midas_output_paths=/u/project/ngarud/Garud_lab/snv_prevalence/scripts/midas_output_paths.txt

OUTDIR=/u/project/ngarud/Garud_lab/snv_prevalence/data

merge_midas.py snps $OUTDIR/snps -i $OUTDIR/midas_output_v1.2.1 -t dir --sample_depth 5 --site_depth 3 --min_samples 1 --max_species 150 --site_prev 0.0 --threads 5 >& $OUTDIR/snps_core_v1.2.1.log

#merge_midas.py snps $OUTDIR/snps_test -i $midas_output_paths -t file --sample_depth 5 --site_depth 3 --min_samples 1 --max_species 150 --site_prev 0.0 --threads 5 >& $OUTDIR/snps_core_v1.2.1.log
