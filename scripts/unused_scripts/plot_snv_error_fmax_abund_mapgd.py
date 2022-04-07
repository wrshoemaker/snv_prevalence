from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import scipy.stats as stats


import diversity_utils
import parse_midas_data
#import calculate_predicted_occupancy
import calculate_predicted_prevalence_mapgd

import prevalence_utils

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import figure_utils



data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])


prevalence_dict = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict.keys())
species_list.sort()

species_list = ['Eubacterium_rectale_56927']

good_species = 'Eubacterium_rectale_56927'
#bad_species = 'Bacteroides_vulgatus_57955'

#clade_types = ['all','largest_clade', 'no_strains']
#pi_types = ['pi_exclude_boundary', 'pi_include_boundary']

clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'


def null_increase_test(relative_error_mean, n_perm=10000):

    relative_error_mean_ratio = relative_error_mean[1:] / relative_error_mean[:-1]
    relative_error_mean_ratio_mean = numpy.mean(relative_error_mean_ratio)

    null_mean_ratio = []
    for i in range(n_perm):

        relative_error_mean_permute = numpy.random.permutation(relative_error_mean)
        relative_error_mean_permute_ratio = relative_error_mean_permute[1:] / relative_error_mean_permute[:-1]
        null_mean_ratio.append(numpy.mean(relative_error_mean_permute_ratio))

    null_mean_ratio = numpy.asarray(null_mean_ratio)

    p_value = sum(null_mean_ratio>relative_error_mean_ratio_mean) / n_perm

    print(relative_error_mean_ratio_mean, p_value)




fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.15)

#ax_f_max = plt.subplot2grid((1, 2), (0,0))
#ax_abundance = plt.subplot2grid((1, 2), (0,1))

ax_f_max = plt.subplot2grid((1, 2), (0,0))
ax_delta_hist = plt.subplot2grid((1, 2), (0,1))

relative_error_ratio_all = []
for species_name in species_list:

    predicted_prevalence = prevalence_dict[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd']
    predicted_prevalence = numpy.asarray(predicted_prevalence)

    observed_prevalence = prevalence_dict[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
    observed_prevalence = numpy.asarray(observed_prevalence)

    f_max = prevalence_dict[species_name][clade_type][pi_type][variant_type]['f_max_mapgd']
    f_max = numpy.asarray(f_max)

    predicted_prevalence_no_zeros = predicted_prevalence[(observed_prevalence>0) & (predicted_prevalence>0)]
    observed_prevalence_no_zeros = observed_prevalence[(observed_prevalence>0) & (predicted_prevalence>0)]
    f_max_no_zeros = f_max[(observed_prevalence>0) & (predicted_prevalence>0)]

    relative_error = numpy.absolute(observed_prevalence_no_zeros - predicted_prevalence_no_zeros) / observed_prevalence_no_zeros

    f_max_no_zeros_to_keep = f_max_no_zeros[f_max_no_zeros<0.5]
    relative_error_to_keep = relative_error[f_max_no_zeros<0.5]


    f_max_no_zeros_to_keep_log10 = numpy.log10(f_max_no_zeros_to_keep)
    # relationship between error and f_max
    hist, bin_edges = numpy.histogram(f_max_no_zeros_to_keep_log10, density=True, bins=20)
    bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]
    #prevalence_predicted_mean = []
    relative_error_to_keep_mean = []
    for i in range(0, len(bin_edges)-1):
        idx_i = (f_max_no_zeros_to_keep_log10 > bin_edges[i]) & (f_max_no_zeros_to_keep_log10 < bin_edges[i+1])
        #prevalence_predicted_mean.append(numpy.mean(numpy.log10(predicted_prevalence[idx_i])))
        relative_error_to_keep_mean.append(numpy.mean(relative_error_to_keep[idx_i]))

    relative_error_to_keep_mean = numpy.asarray(relative_error_to_keep_mean)
    bins_mean = numpy.asarray(bins_mean)

    bins_mean_no_nan = bins_mean[~numpy.isnan(relative_error_to_keep_mean)]
    relative_error_mean_no_nan = relative_error_to_keep_mean[~numpy.isnan(relative_error_to_keep_mean)]


    #relative_error_mean_no_nan

    #null_increase_test(relative_error_mean_no_nan)

    #relative_error_ratio_all.extend(relative_error_mean_no_nan_ratio.tolist())
    relative_error_ratio_all.extend(f_max_no_zeros.tolist())

    #ax_f_max.plot(10**bins_mean_no_nan, relative_error_mean_no_nan, ls ='-', c='k', lw=2, alpha=0.8)

    ax_f_max.scatter(f_max_no_zeros, observed_prevalence_no_zeros, alpha=0.4)





ax_f_max.set_xscale('log', basex=10)
ax_f_max.set_yscale('log', basey=10)




ax_delta_hist.hist(numpy.log10(relative_error_ratio_all), bins=20)


fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%serror_fmax_abund_mapgd.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
