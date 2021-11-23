#!/bin/bash
#$ -N test_HMP_snps
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_HMP_snps_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_HMP_snps_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=60:00:00
#$ -l h_data=16G
#$ -l highp


. /u/local/Modules/default/init/modules.sh
module load python/2.7

export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2

module load samtools
module load bowtie2
module load mapgd



/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/Bacteroides_vulgatus_57955

samtools mpileup -q 20 -Q 25 -B *_sorted.bam \
| mapgd proview -H NZ_GG703852_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered





samtools mpileup -q 20 -Q 25 -B 700175331_sorted.bam 700175086_sorted.bam \
| mapgd proview -H 700175086_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered


samtools mpileup -q 20 -Q 25 -B 700175331_sorted.bam 700173377_sorted.bam \
| mapgd proview -H /u/home/w/wrshoema/project-ngarud/snv_prevalence/scripts/qsub_HMP_MAPGD/test_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered







samtools mpileup -q 5 -Q 5 -B 700175086_sorted.bam \
| mapgd proview -H 700175086_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered.pol


samtools mpileup -q 5 -Q 5 -B 700173377_sorted.bam \
| mapgd proview -H 700173377_sorted.header | mapgd pool -a 20 -o allelefrequency-filtered_700173377_sorted



data_types/file_index.cc:70: no scaffold with name=NZ_GG695952 exists.
./io/map_file.h:776 could not initilize 10File_index

./io/map_file.h:818 could not initilize 5Locus



for d in /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/*/ ; do
    echo "$d"

    while read acc; do
      echo ${acc}
      bam=${d}${acc}_sorted.bam
      header=${d}${acc}_sorted.header
      pol=${d}${acc}_sorted
      if [[ ! -s ${bam} ]] ; then
          echo ${bam} is empty
              continue
      else
        echo ${bam}

        samtools mpileup -q 5 -Q 5 -B ${bam} \
        | mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}

      fi

    done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/HMP1-2_samples_final.txt

done






for d in /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/*/ ; do
    echo "$d"
done




while read acc; do
  echo ${acc}
  bam=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/Bacteroides_vulgatus_57955/${acc}_sorted.bam
  header=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/Bacteroides_vulgatus_57955/${acc}_sorted.header
  pol=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/Bacteroides_vulgatus_57955/${acc}_sorted
  if [[ ! -s ${bam} ]] ; then
      echo ${bam} is empty
          continue
  else
    echo ${bam}

    samtools mpileup -q 5 -Q 5 -B ${bam} \
    | mapgd proview -H ${header} | mapgd pool -a 20 -o ${pol}

  fi

done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/HMP1-2_samples_final.txt
