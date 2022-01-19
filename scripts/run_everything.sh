#!/bin/bash

# activate your conda environment
#source activate microbiome_evolution

# Fig. 1 and S2
python plot_diversity_summary.py

# Fig. S1
python plot_maf_coverage.py

# Figs. S3-S6


python plot_f_mean_vs_prevalence.py
