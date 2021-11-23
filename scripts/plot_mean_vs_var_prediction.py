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


clade_type = 'all'
#clade_type = 'all'
# 'largest_clade'
pi_type = 'pi_include_boundary'
variant_type = '4D'
max_n_occurances = 7

if clade_type == 'all':
    min_prevalence = 0.1
else:
    min_prevalence = 0.1


# sort species names
species_n_sites_list = []
for species_name in species_list:
    observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
    observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)
    n_sites = sum(observed_prevalence_mapgd >= min_prevalence)
    species_n_sites_list.append((species_name, n_sites))

species_n_sites_list.sort(key=lambda tup: tup[1])
species_n_sites_list = species_n_sites_list[::-1]
species_list_to_plot = [s[0] for s in species_n_sites_list]


fig, ax = plt.subplots(figsize=(4,4))

for species_name in species_list_to_plot:

    predicted_f_mean = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_mean_mapgd']
    # predicted_f_var_mapgd
    predicted_f_mean = numpy.asarray(predicted_f_mean)

    predicted_f_var = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_f_mean_mapgd']
    predicted_f_var = numpy.asarray(predicted_f_var)

    observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
    observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)

    if len(predicted_f_var) == 0:
        continue

    predicted_f_mean_to_plot = predicted_f_mean[observed_prevalence_mapgd >= min_prevalence]
    predicted_f_var_to_plot = predicted_f_var[observed_prevalence_mapgd >= min_prevalence]

    ax.scatter(predicted_f_mean_to_plot, predicted_f_var_to_plot, alpha=0.4, s=12, c=species_color_map[species_name])


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.plot([0.001,0.4],[0.001,0.4], lw=3,ls='--',c='k',zorder=1)


fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.22)
fig.savefig("%spredicted_mean_vs_var.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
