#!/bin/bash
#$ -N Alistipes_senegalensis_58364
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/test_mapgd_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=96:00:00
#$ -l h_data=80G
#$ -l highp

#module load python/2.7

. /u/local/Modules/default/init/modules.sh
module use /u/project/ngarud/apps/modulefiles
module load mapgd
module load samtools

species=Alistipes_finegoldii_56071



samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam
| mapgd proview -H seq1.header | mapgd pool -a 22 -o allelefrequency-filtered.pol



# qrsh -l h_rt=3:00:00,h_data=24G



samtools mpileup -q 25 -Q 25 -B 700111439_sorted.bam 700117828c_sorted.bam \
| mapgd proview -H 700111439_sorted.header | mapgd pool -a 22 -o test.pol




samtools mpileup -q 25 -Q 25 -B 700101840_sorted.bam 700109621_sorted.bam \
| mapgd proview -H 700101840_sorted.header | mapgd pool -a 22 -o test.pol



mapgd proview -H 700101840_sorted.header
