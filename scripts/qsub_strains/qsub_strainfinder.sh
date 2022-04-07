#!/bin/bash

#$ -N af_rest
#$ -e /u/project/ngarud/wrshoema/snv_prevalence/strainfinder_input/af_$TASK_ID
#$ -o /u/project/ngarud/wrshoema/snv_prevalence/strainfinder_input/af_$TASK_ID
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=330:00:00
#$ -l h_data=24G
#$ -l highp
#$ -t 1-164


#SAMPLE=$(sed "${SGE_TASK_ID}q;d" "hmp_sample_names.txt")
SAMPLE=$(sed "${SGE_TASK_ID}q;d" "metadata/HMP1-2_sample_ids.txt")
DIR_PATH="/u/project/ngarud/wrshoema/snv_prevalence/data/midas_output_v1.2.1/"

SPECIES_LIST=('Alistipes_finegoldii_56071' 'Alistipes_onderdonkii_55464' 'Alistipes_putredinis_61533' 'Alistipes_shahii_62199' 'Bacteroidales_bacterium_58650' 'Bacteroides_caccae_53434' 'Bacteroides_cellulosilyticus_58046' 'Bacteroides_fragilis_54507' 'Bacteroides_ovatus_58035' 'Bacteroides_stercoris_56735' 'Bacteroides_thetaiotaomicron_56941' 'Bacteroides_uniformis_57318' 'Bacteroides_vulgatus_57955' 'Bacteroides_xylanisolvens_57185' 'Barnesiella_intestinihominis_62208' 'Dialister_invisus_61905' 'Eubacterium_rectale_56927' 'Oscillibacter_sp_60799' 'Parabacteroides_distasonis_56985' 'Parabacteroides_merdae_56972' 'Ruminococcus_bicirculans_59300' 'Ruminococcus_bromii_62047')


#echo ${SAMPLE}
#echo $DIR_PATH
#list of files in form {SPECIES}.snps.gz
#SPECIES_LIST=$(ls $DIR_PATH${SAMPLE}/snps/output/)

#i=0
#while read line
#do
#    SPECIES_LIST[ $i ]="$line"
#    (( i++ ))
#done < <(ls $DIR_PATH${SAMPLE}/snps/output/)

#echo ${SPECIES_LIST}
#echo "______END OF SPECIES LIST_______"

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15 ##--enable-unicode=ucs4
##./configure --enable-unicode=ucs4
####source activate midas
pip uninstall -y scipy
pip install scipy
pip install openopt
pip install FuncDesigner
pip install DerApproximator

START_SPECIES=${SPECIES_LIST[0]}

for index in "${SPECIES_LIST[@]}"
do
    START_SPECIES=$index
    python /u/project/ngarud/wrshoema/snv_prevalence/scripts/get_strainfinder_input.py  ${START_SPECIES}  ${SAMPLE}  /u/project/ngarud/wrshoema/snv_prevalence/strainfinder_input/
    python /u/project/ngarud/wrshoema/snv_prevalence/scripts/run_strainfinder.py  ${START_SPECIES}  ${SAMPLE}
    python /u/project/ngarud/wrshoema/snv_prevalence/scripts/postprocess_strainfinder_output.py  ${START_SPECIES}  ${SAMPLE}
done
