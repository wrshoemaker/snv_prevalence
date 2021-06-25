#!/bin/bash
#$ -N Anaerostipes_hadrus_55206
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=96:00:00
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


species=Anaerostipes_hadrus_55206
#species=$1


dir=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/

#rm /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/*.pol
#rm /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/*.idx

# remove species with no samples
#for d in /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/*/ ; do
#    rmdir ${d}
#done


# do a test run, figure this shit out


#bam=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/Bacteroides_vulgatus_57955/700175331_sorted.bam
#header=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/Bacteroides_vulgatus_57955/700175331_sorted.header
#pol=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/Bacteroides_vulgatus_57955/700175331_sorted


#samtools mpileup -q 5 -Q 5 -B ${bam} \
#| mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}




#for d in /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/*/ ; do
#    echo "$d"

#    while read acc; do
#      echo ${acc}
#      bam=${d}${acc}_sorted.bam
#      header=${d}${acc}_sorted.header
#      pol=${d}${acc}_sorted
#      if [[ ! -s ${bam} ]] ; then
#          echo ${bam} is empty
#              continue
#      else
#        echo ${bam}

#        samtools mpileup -q 5 -Q 5 -B ${bam} \
#        | mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}

#      fi

#    done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/HMP1-2_samples_final.txt

#done




while read acc; do
  echo ${acc}
  bam=${dir}${acc}_sorted.bam
  header=${dir}${acc}_sorted.header
  pol=${dir}${acc}_sorted
  if [[ ! -s ${bam} ]] ; then
      echo ${bam} is empty
          continue

  else
    echo ${bam}

    samtools mpileup -q 5 -Q 5 -B ${bam} \
    | mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}

  fi

done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/samples_species/HMP1-2_samples_final_${species}.txt
