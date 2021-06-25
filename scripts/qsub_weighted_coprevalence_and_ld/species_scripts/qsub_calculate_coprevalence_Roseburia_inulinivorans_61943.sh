#!/bin/bash
#$ -N calculate_coprevalence
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence_output
#$ -l h_data=8G
#$ -l time=144:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13



#python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence_f0.py --ld_moment 1 --species Roseburia_inulinivorans_61943
#python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_coprevalence_f0.py --ld_moment 2 --species Roseburia_inulinivorans_61943


#python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria_f0.py --ld_moment 1 --species Roseburia_inulinivorans_61943
python /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria_f0.py --ld_moment 2 --species Roseburia_inulinivorans_61943
