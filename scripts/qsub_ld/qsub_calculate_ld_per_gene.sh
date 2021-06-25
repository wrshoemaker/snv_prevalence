#!/bin/bash
#$ -N ld_per_gene
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_ld_per_gene_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_ld_per_gene_output
#$ -l h_data=8G
#$ -l time=36:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
#module load python/3.6.1
module load python/2.7.13


python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria_per_gene.py
