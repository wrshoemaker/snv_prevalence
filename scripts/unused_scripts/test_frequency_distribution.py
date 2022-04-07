from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path
import random

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
import scipy.stats as stats
import scipy.special as special

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import calculate_predicted_prevalence_mapgd

import plot_utils

species_name = 'Bacteroides_ovatus_58035'
variant_type = '4D'

def get_f_no_zeros_dict(species_name, clade_type='all'):

    intermediate_filename_template = config.data_directory+"frequency_dict_mapgd/%s_%s.dat"
    intermediate_filename = intermediate_filename_template % (species_name, clade_type)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b


f_no_zeros_mapgd_dict = get_f_no_zeros_dict(species_name)
#sf_no_zeros_mapgd = f_no_zeros_mapgd_dict[variant_type]

sites = list(f_no_zeros_mapgd_dict[variant_type].keys())

# remove sites with fixed.....

random.shuffle(sites)

sites_range = range(len(sites))
sites_to_plot_chunked = [sites[x:x+50] for x in range(0, len(sites_range), 50)]


#print(sites_to_plot_chunked)

fig, ax = plt.subplots(figsize=(4,4))



for sites_to_plot in sites_to_plot_chunked:

    all_freqs = []

    for s in sites_to_plot:

        f_no_zeros_mapgd = f_no_zeros_mapgd_dict[variant_type][s]['frequencies']
        all_freqs.extend(f_no_zeros_mapgd)

    f_no_zeros_mapgd = numpy.asarray(all_freqs)
    f_no_zeros_mapgd = f_no_zeros_mapgd[ (f_no_zeros_mapgd>0) & (f_no_zeros_mapgd<1)]
    f_no_zeros_mapgd_log10 = numpy.log10(f_no_zeros_mapgd)
    #f_no_zeros_mapgd_log10_rescaled = (f_no_zeros_mapgd_log10 - numpy.mean(f_no_zeros_mapgd_log10)) / numpy.std(f_no_zeros_mapgd_log10)
    hist_f, bin_edges_f = numpy.histogram(f_no_zeros_mapgd_log10, density=True, bins=20)
    bins_mean_f = [0.5 * (bin_edges_f[i] + bin_edges_f[i+1]) for i in range(0, len(bin_edges_f)-1 )]
    bins_mean_f = numpy.asarray(bins_mean_f)
    hist_f_to_plot = hist_f[hist_f>0]
    bins_mean_f_to_plot = bins_mean_f[hist_f>0]

    ax.plot(10**bins_mean_f_to_plot, hist_f_to_plot, ls='-', c='dodgerblue', alpha=0.05, lw=2)


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel('Allele frequencies', fontsize=12)
ax.set_ylabel('Probability density', fontsize=12)

fig.tight_layout()
fig.subplots_adjust(wspace=0.34, hspace=0.28)
fig.savefig("%stest_frequency_distribution.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()

#['frequencies']
#print(f_no_zeros_mapgd[('NC_009614', 3520674L)]['frequencies'])
