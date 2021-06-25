#!/bin/bash
#$ -N run_ld_all_taxa
#$ -e /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/run_ld_all_taxa_error
#$ -o /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/run_ld_all_taxa_output
#$ -l h_data=8G
#$ -l time=6:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.13


python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.0 --high 0.05

python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.0 --high 0.1

python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.0 --high 0.2

python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.0 --high 0.4

python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.0 --high 0.6

python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.0 --high 0.8

#python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.0 --high 0.8


#python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.05 --high 1

#python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.1 --high 1

#python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.3 --high 1

#python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.5 --high 1

#python2 /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py --condition_on_freq --low 0.7 --high 1


#python2 /u/home/w/wrshoema/project-ngarud/negative_selection_microbiome/scripts/calculate_linkage_disequilibria.py



#rsync -avR wrshoema@hoffman2.idre.ucla.edu:/u/home/w/wrshoema/project-ngarud/negative_selection_microbiome/data/snps/\*/coverage_distribution.txt.bz2 .
