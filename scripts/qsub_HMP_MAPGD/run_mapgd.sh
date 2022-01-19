#!/bin/bash
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=96:00:00
#$ -l h_data=80G
#$ -l highp

#module load python/2.7




#####  $ -N test_mapgd

#export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
#export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
#export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2


. /u/local/Modules/default/init/modules.sh
module use /u/project/ngarud/apps/modulefiles
module load mapgd
module load samtools


#species=Bacteroides_vulgatus_57955
species=$1


dir=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output_bam/${species}/

#rm /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/*.pol
#rm /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/*.idx



#while read acc; do
#  echo ${acc}
#  bam=${dir}${acc}_sorted.bam
#  header=${dir}${acc}_sorted.header
#  pol=${dir}${acc}_sorted
#  if [[ ! -s ${bam} ]] ; then
#      echo ${bam} is empty
#          continue

#  else
#    echo ${bam}

    #samtools mpileup -q 5 -Q 5 -B ${bam} \
    #| mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}

#  fi

#done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/samples_species/HMP1-2_samples_final_${species}.txt




    #samtools mpileup -q 5 -Q 5 -B ${bam} \
    #| mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}



# add option to ignore files smaller than






#while read acc; do
#  echo ${acc}
#  bam=${dir}${acc}_sorted.bam
#  header=${dir}${acc}_sorted.header
#  pol=${dir}${acc}_sorted
#  if [ ! -s ${bam} ] || [ stat -c%s ${bam} -le 50000 ]; then
#      echo ${acc} is empty
#          continue

#  else
#    echo ${acc}

#    #samtools mpileup -q 5 -Q 5 -B ${bam} \
#    #| mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}

#  fi

#done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/samples_species/HMP1-2_samples_final.txt


while read acc; do
  bam=${dir}${acc}_sorted.bam
  header=${dir}${acc}_sorted.header
  pol=${dir}${acc}_sorted
  if [[ ! -s ${bam} ]]; then
      echo ${acc} is empty
          continue

  else
    filesize=$(stat -c%s "$bam")
    if [ $filesize -le 1000000 ]; then
        echo ${acc} is too small
            continue
    else
        echo ${acc}
        samtools mpileup -q 5 -Q 5 -B ${bam} \
        | mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}

    fi

  fi

done </u/project/ngarud/wrshoema/snv_prevalence/scripts/qsub_HMP_MAPGD/samples_species/HMP1-2_samples_final.txt
