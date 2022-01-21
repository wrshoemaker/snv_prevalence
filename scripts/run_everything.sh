#!/bin/bash

# activate your conda environment
#source activate microbiome_evolution

# Fig. 1 and S2
python plot_diversity_summary.py

# Fig. 2
python plot_prevalence_prediction_summary.py

# Fig. 3
python plot_prevalence_vs_error_mapgd_slm.py


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


python plot_f_mean_vs_prevalence.py


python plot_prevalence_vs_error_mapgd_slm.py
