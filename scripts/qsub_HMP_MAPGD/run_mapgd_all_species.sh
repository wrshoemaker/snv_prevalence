#!/bin/bash



. /u/local/Modules/default/init/modules.sh
module load python/2.7



for d in /u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/*/ ; do
  species="$(echo "$d" | rev | cut -d '/' -f 2 | rev)"
  echo ${species}
  qsub -N ${species} /u/home/w/wrshoema/project-ngarud/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/run_mapgd.sh ${species}
done


#species=Anaerostipes_hadrus_55206

#qsub -N ${species} /u/home/w/wrshoema/project-ngarud/negative_selection_microbiome/scripts/qsub_HMP_MAPGD/run_mapgd.sh ${species}
