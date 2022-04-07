#!/bin/bash

# activate your conda environment
#source activate microbiome_evolution

# run MIDAS
#sh ./qsub_midas/qsub_script_HMP_MIDAS_species_snps_cnvs.sh

# merge MIDAS species
#sh ./qsub_midas/qsub_script_HMP_MIDAS_species_merge.sh

# merge MIDAS genes
#sh ./qsub_midas/qsub_script_HMP_MIDAS_CNVs_merge.sh

# merge MIDAS SNPs
#sh ./qsub_midas/qsub_script_HMP_MIDAS_SNPs_merge.sh

# postproces midas data
#python postprocess_all_midas_data_serial.py

# split all BAM files into separate species
#sh ./qsub_HMP_MAPGD/split_all_bam_files.sh

# run MAPGD on each host for each species
#sh ./qsub_HMP_MAPGD/run_mapgd_all_species.sh

# make dictionary of alleles for all species
#sh ./qsub_prevalence/run_make_alleles_dict_all_species.sh

# predict prevalence for all species
#sh ./qsub_prevalence/run_prevalence_prediction_all_species.sh

# run StrainFinder
#sh ./qsub_strains/qsub_strainfinder.sh



# Fig. 1 and S2
python plot_diversity_summary.py

# Fig. 2
python plot_prevalence_prediction_example.py

# Fig. 3
python plot_prevalence_prediction_summary.py

# Fig. 4
#python plot_prevalence_vs_error_mapgd_slm.py
python plot_slm_error_vs_strainfinder

# Supplement
# Fig. S1
python plot_maf_coverage.py

# Figs. S2, S3
python plot_f_mean_vs_beta.py


# Figs. S4-S7

python plot_predicted_observed_f_mean.py
python plot_predicted_observed_f_var.py

python plot_error_evo_vs_slm.py

python plot_prevalence_error_mapgd.py

python python plot_predicted_observed_prevalence_mapgd.py

python plot_f_mean_vs_prevalence.py

python plot_f_max_vs_prevalence.py


#python plot_prevalence_vs_error_mapgd_slm.py
