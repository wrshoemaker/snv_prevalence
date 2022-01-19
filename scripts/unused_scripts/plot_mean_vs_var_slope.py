from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
import scipy.stats as stats

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import calculate_predicted_prevalence_mapgd
import calculate_predicted_prevalence


import plot_utils

species_color_map, ordered_species_list = plot_utils.get_species_color_map()


prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()


#fig, ax = plt.subplots(figsize=(4,4))

clade_type = 'no_strains'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'
max_n_occurances = 7


if clade_type == 'all':
    min_prevalence = 0.2
else:
    min_prevalence = 0.1


means_all = []
variances_all = []

predicted_means_all = []
predicted_variances_all = []

for species_name in species_list:

    f_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_no_zeros_mapgd']
    f_mapgd = numpy.asarray(f_mapgd)

    n_non_zero_f = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
    n_non_zero_f = numpy.asarray(n_non_zero_f)

    f_mean_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_mean_mapgd']
    f_mean_mapgd = numpy.asarray(f_mean_mapgd)

    f_var_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_var_mapgd']
    f_var_mapgd = numpy.asarray(f_var_mapgd)

    f_max_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_max_mapgd']
    f_max_mapgd = numpy.asarray(f_max_mapgd)

    predicted_f_mean = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_f_mean_mapgd']
    predicted_f_mean = numpy.asarray(predicted_f_mean)

    predicted_f_var = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_f_var_mapgd']
    predicted_f_var = numpy.asarray(predicted_f_var)

    observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
    observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)

    if len(f_mapgd) == 0:
        continue

    f_mean_mapgd_to_plot = f_mean_mapgd[observed_prevalence_mapgd >= min_prevalence]
    f_var_mapgd_to_plot = f_var_mapgd[observed_prevalence_mapgd >= min_prevalence]

    predicted_f_mean_mapgd_to_plot = predicted_f_mean[observed_prevalence_mapgd >= min_prevalence]
    predicted_f_var_mapgd_to_plot = predicted_f_var[observed_prevalence_mapgd >= min_prevalence]


    means_all.extend(f_mean_mapgd_to_plot.tolist())
    variances_all.extend(f_var_mapgd_to_plot.tolist())


    predicted_means_all.extend(predicted_f_mean_mapgd_to_plot.tolist())
    predicted_variances_all.extend(predicted_f_var_mapgd_to_plot.tolist())





#slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(means_all), numpy.log10(variances_all))
#slope_predicted, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(predicted_means_all), numpy.log10(predicted_variances_all))


fig, ax = plt.subplots(figsize=(4,4))
fig.subplots_adjust(bottom= 0.15)

ax.scatter(predicted_means_all, predicted_variances_all, s=15, alpha=0.8)

ax.set_xlim(0.00001, 1)
ax.set_ylim(0.00001, 1)

ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)




fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.25)
fig.savefig("%spredicted_mean_vs_var.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
