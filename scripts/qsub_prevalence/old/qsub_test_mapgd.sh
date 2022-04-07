#!/bin/bash
#$ -N MIDAS_1
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/postproc_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/postproc_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=120:00:00
#$ -l h_data=34G
#$ -l highp

####ttt##$-tc 100 # Throttle to max 100 tasks at a time
########$-t 1-158
#######$-ltime=23:00:00


. /u/local/Modules/default/init/modules.sh


module load samtools
module load mapgd
module load bowtie2



# paried reads only
#samtools fastq -f 0x1 -F 0xC NZ_GG703852.bam -1 NZ_GG703852_R1.fastq.gz -2 NZ_GG703852_R2.fastq.gz

samtools fastq -f 0x02 NZ_GG703852.bam -1 NZ_GG703852_R1.fastq.gz -2 NZ_GG703852_R2.fastq.gz



bowtie2-build genome.fna Prevotella_copri_61740


#bowtie2 -x ./Prevotella_copri_61740 \
#  -1 NZ_GG703852_R1.fastq.gz -2 NZ_GG703852_R2.fastq.gz -S test.sam


# thank god
bowtie2 -x ./Prevotella_copri_61740 \
  -1 NZ_GG703852_R1.fastq.gz -2 NZ_GG703852_R2.fastq.gz \
  | samtools view -bS - > NZ_GG703852_test.bam


samtools sort -T /tmp/aln.sorted -o NZ_GG703852_test_sorted.bam NZ_GG703852_test.bam

samtools view -H NZ_GG703852_test_sorted.bam > seq1.header


samtools mpileup -q 25 -Q 25 -B NZ_GG703852_test_sorted.bam \
| mapgd proview -H seq1.header | mapgd pool -a 22 -o allelefrequency-filtered





# try with original bam

samtools sort -T /tmp/aln.sorted -o NZ_GG703852_sorted.bam NZ_GG703852.bam
samtools view -H NZ_GG703852_sorted.bam > NZ_GG703852_sorted.header


samtools mpileup -q 20 -Q 25 -B NZ_GG703852_sorted.bam \
| mapgd proview -H NZ_GG703852_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered




samtools mpileup -q 5 -Q 5 -B 700175086_sorted.bam 700175181_sorted.bam \
| mapgd proview -H 700175181_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered.pol




samtools mpileup -q 5 -Q 5 -B 700173377_sorted.bam 700175086_sorted.bam \
| mapgd proview -H /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/test_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered.pol



samtools mpileup -q 5 -Q 5 -B 700173377_sorted.bam 700175086_sorted.bam \
| mapgd proview -H 700173377_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered.pol




samtools mpileup -q 5 -Q 5 -B 700175331_sorted.bam \
| mapgd proview -H 700175331_sorted.header | mapgd pool -a 22 -o allelefrequency-filtered_700175331_sorted






#/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/rep_genomes/Bacteroides_vulgatus_57955/genome.fna.fai
# remove temp files






#while read acc; do
#  echo ${acc}
#  rm -r /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output/${acc}/snps/temp/
#done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/HMP1-2_samples_final.txt
