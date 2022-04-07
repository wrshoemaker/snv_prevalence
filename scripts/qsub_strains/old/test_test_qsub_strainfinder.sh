#!/bin/bash

#$ -N rest_af
#$ -e /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_input/af$TASK_ID
#$ -o /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_input/af$TASK_ID
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=72:00:00
#$ -l h_data=12G
#$ -l highp
#$ -t 1-164


#############ONLY CHANGED TEMPORARILY TO COMPARE W RICKYS DATA
#SAMPLE=$(sed "${SGE_TASK_ID}q;d" "hmp_sample_names.txt")
SAMPLE=$(sed "${SGE_TASK_ID}q;d" "sample_names/sample_list.txt")
DIR_PATH="/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/midas_output_v1.2.1/"
#DIR_PATH="/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output/"

echo ${SAMPLE}
#echo $DIR_PATH
#list of files in form {SPECIES}.snps.gz
#SPECIES_LIST=$(ls $DIR_PATH${SAMPLE}/snps/output/)

i=0
while read line
do
    SPECIES_LIST[ $i ]="$line"
    (( i++ ))
done < <(ls $DIR_PATH${SAMPLE}/snps/output/)

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

    python /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/test_af_get_strainfinder_input.py  ${START_SPECIES}  ${SAMPLE}  /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_input/

##    python /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/test_ricky_get_strainfinder_input.py  ${START_SPECIES}  ${SAMPLE}  /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_input/

    python /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/test_ricky_strainfinder.py  ${START_SPECIES}  ${SAMPLE}

    python /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/test_ricky_postprocess_strainfinder_output.py  ${START_SPECIES}  ${SAMPLE}

done
