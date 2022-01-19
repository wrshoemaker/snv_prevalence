#!/bin/bash
#$ -N get_coverage
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/get_coverage_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/get_coverage_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=48:00:00
#$ -l h_data=24G
#$ -l highp



module unload python
module load anaconda/python2-4.2
source activate midas
module load samtools


#qrsh -l h_rt=3:00:00,h_data=24G

species=Bacteroides_ovatus_58035
#sample=700101916

#ref_genome=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/rep_genomes/${species}/genome.fna



#out=/u/home/w/wrshoema/project-ngarud/snv_prevalence/data/bam_coverage/${species}/${sample}.txt

#samtools mpileup -q 5 -Q 5 -B ${bam} | samtools depth /dev/stdin > test.txt
#samtools mpileup -q 5 -Q 5 -B ${bam} > ${out}


declare -a all_species=('Alistipes_finegoldii_56071' 'Alistipes_onderdonkii_55464' 'Alistipes_putredinis_61533' 'Alistipes_shahii_62199' 'Bacteroidales_bacterium_58650' 'Bacteroides_caccae_53434' 'Bacteroides_cellulosilyticus_58046' 'Bacteroides_fragilis_54507' 'Bacteroides_ovatus_58035' 'Bacteroides_stercoris_56735' 'Bacteroides_thetaiotaomicron_56941' 'Bacteroides_uniformis_57318' 'Bacteroides_vulgatus_57955' 'Bacteroides_xylanisolvens_57185' 'Barnesiella_intestinihominis_62208' 'Dialister_invisus_61905' 'Eubacterium_rectale_56927' 'Oscillibacter_sp_60799' 'Parabacteroides_distasonis_56985' 'Parabacteroides_merdae_56972' 'Ruminococcus_bicirculans_59300' 'Ruminococcus_bromii_62047')


#for species in "${all_species[@]}"
#do
mkdir -p /u/home/w/wrshoema/project-ngarud/snv_prevalence/data/bam_coverage/${species}
for filename in /u/home/w/wrshoema/project-ngarud/snv_prevalence/data/midas_output_bam/${species}/*_sorted.bam; do
    echo $filename
    sample="$(echo "$filename" | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)"
    #echo $sample
    bam=/u/project/ngarud/wrshoema/snv_prevalence/data/midas_output_bam/${species}/${sample}.bam
    out=/u/project/ngarud/wrshoema/snv_prevalence/data/bam_coverage/${species}/${sample}.gz
    # same flag settings used in MAPGD
    samtools mpileup -q 5 -Q 5 -B ${bam} | cut -f 1-2,4 | gzip > ${out}
done
#done








#samtools depth deduped_MA605.bam > .coverage



#samtools mpileup -f ${ref_genome} -r ${bam} | cut -f 5 | tr '[a-z]' '[A-Z]' | fold -w 1 | sort | uniq -c > ${out}
