from __future__ import division
import sample_utils
import config
import parse_midas_data
import os.path
import os
import pylab
import sys
import numpy
import gzip
import pickle
import bz2

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
from math import log10,ceil,fabs,isnan
from numpy.random import randint, choice, multinomial

import midas_db_utils

import scipy.stats as stats
from scipy.stats import gamma, gaussian_kde

import parse_HMP_data
import figure_utils
import calculate_predicted_occupancy

#good_species_list = [good_species_list[3]]

import matplotlib.pyplot as plt
import matplotlib.cm as cm


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()


f_max_all = []

for species_name in prevalence_dict.keys():

    f_max = prevalence_dict[species_name]['all']['4D']['f_max']['to_plot']

    f_max_all.extend(f_max)

f_max_all = numpy.asarray(f_max_all)

f_max_all_log10 = numpy.log10(f_max_all)

f_max_all_log10_mean = numpy.mean(f_max_all_log10)
f_max_all_log10_std = numpy.std(f_max_all_log10)



fig, ax = plt.subplots(figsize=(4,4))

for species_name in prevalence_dict.keys():

    f_max = prevalence_dict[species_name]['all']['4D']['f_max']['to_plot']
    f_max = numpy.asarray(f_max)
    f_max_log10 = numpy.log10(f_max)

    f_max_log10_rescaled = (f_max_log10 - f_max_all_log10_mean) / f_max_all_log10_std

    hist, bin_edges = numpy.histogram(f_max_log10_rescaled, density=True, bins=20)

    bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]

    ax.scatter(bins_mean, hist, alpha=0.8)


ax.set_yscale('log', basey=10)

ax.set_xlabel('Rescaled log ' + r'$f_{max}$', fontsize=14)
ax.set_ylabel('Probability density', fontsize=14)




fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%sf_max_distribution.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
