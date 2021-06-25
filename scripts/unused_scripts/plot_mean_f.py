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

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_occupancy

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])
clade_types = ['all','largest_clade']

prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()


species_list = list(prevalence_dict.keys())

def plot_f_mean():

    f_mean_all = []


    for species_name in species_list:

        f_mean = prevalence_dict[species_name]['all']['4D']['f_mean']['to_plot']


        f_mean_all.extend(f_mean)


    f_mean_all = numpy.log10(f_mean_all)

    mean_f_mean = numpy.mean(f_mean_all)
    std_f_mean = numpy.std(f_mean_all)


    fig, ax = plt.subplots(figsize=(4,4))

    for species_name in species_list:

        f_mean = prevalence_dict[species_name]['all']['4D']['f_mean']['to_plot']

        f_mean = numpy.log10(f_mean)

        f_mean_rescaled = (f_mean  - mean_f_mean) / std_f_mean


        hist, bin_edges = numpy.histogram(f_mean_rescaled, density=True, bins=20)

        bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]

        ax.scatter(bins_mean, hist, alpha=0.8)


    #ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

    ax.set_ylim([0.0001, 1])


    fig.tight_layout()
    #fig.subplots_adjust(hspace=0.2)
    fig.savefig("%sf_mean_hist.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



fig, ax = plt.subplots(figsize=(4,4))

species_list = ['Bacteroides_vulgatus_57955']


for species_name in species_list:


    f_max = prevalence_dict[species_name]['all']['4D']['f_max']['to_plot']

    observed_prevalence = prevalence_dict[species_name]['all']['4D']['observed_prevalence_to_plot']

    ax.scatter(f_max, observed_prevalence, alpha=0.4)

    f_max_range = numpy.linspace(0, 1, num=1000)

    predicted_prevalence = 1 - (1/(1 + 10*f_max_range))**0.3
    predicted_prevalence_2 = 1 - (1/(1 + 100*f_max_range))**0.02

    ax.plot(f_max_range, predicted_prevalence, c='k')

    ax.plot(f_max_range, predicted_prevalence_2, c='k')


    f_max_line = prevalence_dict[species_name]['all']['4D']['f_max_line']
    predicted_prevalence_line = prevalence_dict[species_name]['all']['4D']['predicted_prevalence_line']

    f_max_line = numpy.asarray(f_max_line)
    predicted_prevalence_line = numpy.asarray(predicted_prevalence_line)


    f_max_line = f_max_line[numpy.logical_not(numpy.isnan(predicted_prevalence_line))]
    predicted_prevalence_line = predicted_prevalence_line[numpy.logical_not(numpy.isnan(predicted_prevalence_line))]

    ax.plot(10**f_max_line, 10**predicted_prevalence_line, c='k')


    ax.set_title(figure_utils.get_pretty_species_name(species_name), fontsize=12, fontweight='bold', color='k' )



    #f_mean = prevalence_dict[species_name]['all']['4D']['f_max']['to_plot']

ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

#ax.set_ylim([0.1, 1])

fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%sf_max_vs_prevalence_example.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
