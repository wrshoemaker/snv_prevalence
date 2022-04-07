#!/bin/bash



. /u/local/Modules/default/init/modules.sh
module load python/2.7

export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2

module load samtools
module load bowtie2
module load mapgd


# get reference genomes ready
for d in /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/rep_genomes/*/ ; do
    species_genome=${d}genome.fna
    species_genome_gz=${d}genome.fna.gz
    species_genome_fai=${d}genome.fna.fai
    species_genome_bed=${d}genome.bed
    gzip -dc ${species_genome_gz} >> ${species_genome}
    species="$(echo "$d" | rev | cut -d '/' -f 2 | rev)"
    echo $species
    bowtie2-build ${species_genome} ${d}/${species}
    samtools faidx ${species_genome} -o ${species_genome_fai}
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${species_genome_fai} > ${species_genome_bed}
    species_dir=/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/${species}/
    mkdir -p $species_dir
done



# Done
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

while read acc; do
  echo ${acc}
  qsub /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/test_qsub_HMP1-2_MIDAS_snps.sh ${acc}
done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/HMP1-2_samples_final.txt


#</u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/HMP1-2_samples_final_10.txt
#while read p; do
#  echo "$p"
#done </u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/HMP1-2_samples_final.txt
