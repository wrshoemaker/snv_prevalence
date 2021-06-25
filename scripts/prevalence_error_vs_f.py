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

clade_type_label_dict = {'all': 'All hosts', 'largest_clade': 'Largest clade'}


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict.keys())


fig = plt.figure(figsize = (8, 4)) #

ax_f_mean = plt.subplot2grid((1, 2), (0, 0))
ax_f_max = plt.subplot2grid((1, 2), (0, 1))

ax_list = [ax_f_mean, ax_f_max]
estimates = ['f_mean', 'f_max']

relative_error = prevalence_dict[species_list[0]]['all']['relative_error_to_plot']
relative_error = numpy.asarray(relative_error)

for estimate_idx, estimate in enumerate(estimates):

    ax_estimate = ax_list[estimate_idx]

    estimate_array = prevalence_dict[species_list[0]]['all'][estimate]['to_plot']
    estimate_array = numpy.asarray(estimate_array)

    estimate_slope = prevalence_dict[species_list[0]]['all'][estimate]['slope']
    estimate_intercept = prevalence_dict[species_list[0]]['all'][estimate]['intercept']

    all_ = numpy.concatenate([estimate_array, relative_error])

    # Calculate the point density
    xy = numpy.vstack([estimate_array, relative_error])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = estimate_array[idx], relative_error[idx], z[idx]

    ax_estimate.scatter(x, y, c=z, cmap="Blues", s=30, alpha=0.5, edgecolor='', zorder=1)


    x_log10_range =  numpy.linspace(min(numpy.log10(x)) , max(numpy.log10(x)), 10000)
    y_log10_fit_range = 10 ** (estimate_slope*x_log10_range + estimate_intercept)

    ax_estimate.plot(10**x_log10_range, y_log10_fit_range, c='k', lw=2.5, linestyle='--', zorder=2)

    ax_estimate.set_xscale('log', basex=10)
    ax_estimate.set_yscale('log', basey=10)

    ax_estimate.set_xlim([min(x)*0.8,max(x)*1.1])
    ax_estimate.set_ylim([min(y)*100,max(y)*1.1])



fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%sf_vs_rel_error.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
