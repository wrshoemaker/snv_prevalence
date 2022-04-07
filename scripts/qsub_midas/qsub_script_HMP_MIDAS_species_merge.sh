#!/bin/bash
#$ -N species_merge
#$ -e /u/project/ngarud/Garud_lab/snv_prevalence/scripts/species_error
#$ -o /u/project/ngarud/Garud_lab/snv_prevalence/scripts/species_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=72:00:00
#$ -l h_data=34G
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

# -t = directory containing all samples
#merge_midas.py species $OUTDIR/species -i $OUTDIR/midas_output_v1.2.1 -t dir  >& $OUTDIR/species.log


merge_midas.py species $OUTDIR/species_test -i $midas_output_paths -t file  >& $OUTDIR/species_test.log
