#!/bin/bash
#$ -N run_pairwise_all_taxa
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/run_pairwise_all_taxa_output
#$ -l h_data=8G
#$ -l time=1:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13


python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_distance_vs_gene_overlap.py
