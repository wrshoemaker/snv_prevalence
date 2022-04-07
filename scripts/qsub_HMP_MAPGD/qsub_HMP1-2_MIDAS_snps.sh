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


#export MIDAS_DB=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_db_v1.2

acc=$1




#module load mapgd

#readarray accs < /u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/fix.txt
#accs=(null ${accs[@]}) # zero to one start index
#acc=${accs[$SGE_TASK_ID]}
#echo $acc

#readarray accs < /u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/HMP1-2_samples_final.txt
#OUTDIR=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output/${acc}

#acc=700014954



# get reference genomes ready
#for d in /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/rep_genomes/*/ ; do
#    species_genome=${d}genome.fna
#    species_genome_gz=${d}genome.fna.gz
#    species_genome_fai=${d}genome.fna.fai
#    species_genome_bed=${d}genome.bed
#    gzip -dc ${species_genome_gz} >> ${species_genome}
#    species="$(echo "$d" | rev | cut -d '/' -f 2 | rev)"
#    echo $species
#    bowtie2-build ${species_genome} ${d}/${species}
#    samtools faidx ${species_genome} -o ${species_genome_fai}
#    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${species_genome_fai} > ${species_genome_bed}
#    species_dir=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/
#    mkdir -p ${species_dir}
#done





OUTDIR=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output/${acc}
# Assume OUTDIR already made during species step
fastq1=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/fastq_files/${acc}_1.fastq.gz
fastq2=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/fastq_files/${acc}_2.fastq.gz
species_union=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/species_unions/${acc}_species_union.txt
#species_profile=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output/${acc}/species/species_profile.tx

run_midas.py snps $OUTDIR -1 $fastq1 -2 $fastq2
# --extra_species_file $species_union

# --extra_species_file is causing the problem here
#Error message:
#bowtie2-build: error while loading shared libraries: libtbbmalloc_proxy.so.2: cannot open shared object file: No such file or directory


bam=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output/${acc}/snps/temp/genomes.bam
# index bam
samtools index ${bam}


# create folder for bam output
while read species; do
  echo "$species"
  species_genome_bed=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/rep_genomes/${species}/genome.bed
  species_dir=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output_bam/${species}/
  bowtie2_build=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/rep_genomes/${species}/${species}/

  samtools view -L ${species_genome_bed} -o ${species_dir}${acc}.bam ${bam}
  # delete bam file if its zero

  samtools sort -T /tmp/aln.sorted -o ${species_dir}${acc}_sorted.bam ${species_dir}${acc}.bam
  samtools view -H ${species_dir}${acc}_sorted.bam > ${species_dir}${acc}_sorted.header

done </u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output/${acc}/snps/species.txt



#species_genome=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/rep_genomes/Prevotella_copri_61740/genome.fna
#species_genome_gz=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/rep_genomes/Prevotella_copri_61740/genome.fna.gz

# and then use these sorted bam files for mapgd
