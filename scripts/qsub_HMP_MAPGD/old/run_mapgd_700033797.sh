#!/bin/bash
#$ -N test_mapgd
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=24:00:00
#$ -l h_data=80G
#$ -l highp

#module load python/2.7

#export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
#export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
#export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2


. /u/local/Modules/default/init/modules.sh
module use /u/project/ngarud/apps/modulefiles
module load mapgd
module load samtools



species=Bacteroides_vulgatus_57955
#species=$1


dir=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/

#rm /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/*.pol
#rm /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/*.idx


bam=${dir}700033797_sorted.bam
header=${dir}700033797_sorted.header
pol=${dir}700033797_sorted_test


samtools mpileup -q 25 -Q 25 -B ${bam} \
| mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}
