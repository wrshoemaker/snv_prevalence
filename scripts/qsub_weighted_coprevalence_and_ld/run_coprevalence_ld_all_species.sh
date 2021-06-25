#!/bin/bash



. /u/local/Modules/default/init/modules.sh
#module load python/2.7




while read species; do
  echo "$species"
  qsub -N ${species} /u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_weighted_coprevalence_and_ld/qsub_calculate_coprevalence_and_ld.sh ${species}
done </u/project/ngarud/wrshoema/negative_selection_microbiome/scripts/qsub_weighted_coprevalence_and_ld/good_species_list.txt
