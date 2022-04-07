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

import parse_midas_data

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_occupancy

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])



prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict.keys())
species_list.sort()


fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.15)

ax_f_max = plt.subplot2grid((1, 2), (0,0))
ax_f_mean = plt.subplot2grid((1, 2), (0,1))


for species_name in species_list:

    f_max_line = prevalence_dict[species_name]['all']['4D']['f_max_line']
    f_max_line = numpy.asarray(f_max_line)

    f_mean_line = prevalence_dict[species_name]['all']['4D']['f_mean_line']
    f_mean_line = numpy.asarray(f_mean_line)

    f_max_relative_error_line = prevalence_dict[species_name]['all']['4D']['f_max_vs_relative_error_line']
    f_mean_relative_error_line = prevalence_dict[species_name]['all']['4D']['f_mean_vs_relative_error_line']

    ax_f_max.plot(10**f_max_line, f_max_relative_error_line, ls ='-', c='k', alpha=0.4)
    ax_f_mean.plot(10**f_mean_line, f_mean_relative_error_line, ls ='-', c='k', alpha=0.4)




ax_f_max.set_xscale('log', basex=10)
ax_f_mean.set_xscale('log', basex=10)


ax_f_max.set_xlabel('Maximum frequency\nacross hosts, ' + r'$f_{max}$', fontsize=14)

ax_f_mean.set_xlabel('Mean frequency\nacross hosts, ' + r'$\left \langle f \right \rangle$', fontsize=14)

ax_f_max.set_ylabel('SNV prevalence relative error', fontsize=14)
ax_f_mean.set_ylabel('SNV prevalence relative error', fontsize=14)

fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
fig.savefig("%sfmax_vs_error.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
