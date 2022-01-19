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



### -N ${acc}
### -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/qsub_out/${acc}_error
### -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/qsub_out/${acc}_output

. /u/local/Modules/default/init/modules.sh
#module load python/2.7


module unload python
module load anaconda/python2-4.2

source activate midas

module load samtools
#module load bowtie2

export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2


acc=700013588

species=Bacteroides_ovatus_58035

species_genome_bed=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/rep_genomes/${species}/genome.bed
species_dir=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output_bam/${species}/
bowtie2_build=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/rep_genomes/${species}/${species}/


bam=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/data/midas_output/${acc}/snps/temp/genomes.bam
#echo $acc
samtools view -L ${species_genome_bed} -o ${species_dir}${acc}.bam ${bam}
# delete bam file if its zero
samtools sort -T /tmp/aln.sorted -o ${species_dir}${acc}_sorted.bam ${species_dir}${acc}.bam
samtools view -H ${species_dir}${acc}_sorted.bam > ${species_dir}${acc}_sorted.header


while read acc; do
  bam=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output/${acc}/snps/temp/genomes.bam
  #echo $acc
  samtools view -L ${species_genome_bed} -o ${species_dir}${acc}.bam ${bam}
  # delete bam file if its zero
  samtools sort -T /tmp/aln.sorted -o ${species_dir}${acc}_sorted.bam ${species_dir}${acc}.bam
  samtools view -H ${species_dir}${acc}_sorted.bam > ${species_dir}${acc}_sorted.header
done </u/project/ngarud/wrshoema/snv_prevalence/scripts/qsub_HMP_MAPGD/HMP1-2_samples_final.txt
