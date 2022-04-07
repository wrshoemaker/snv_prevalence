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
#clade_type = 'largest_clade'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'
#max_n_occurances = 7
n_sample = 1000
min_n_non_zero_f = 3


fig = plt.figure(figsize = (8.5, 15)) #
#fig.subplots_adjust(bottom= 0.1)

ax_mean = plt.subplot2grid((4, 2), (0, 0), colspan=1)
ax_mean_error = plt.subplot2grid((4, 2), (0, 1), colspan=1)

ax_var = plt.subplot2grid((4, 2), (1, 0), colspan=1)
ax_var_error = plt.subplot2grid((4, 2), (1, 1), colspan=1)

ax_prevalence = plt.subplot2grid((4, 2), (2, 0), colspan=1)
ax_prevalence_error = plt.subplot2grid((4, 2), (2, 1), colspan=1)

ax_f_max_vs_prevalence = plt.subplot2grid((4, 2), (3, 0), colspan=1)
ax_f_max_vs_prevalence_error = plt.subplot2grid((4, 2), (3, 1), colspan=1)



if clade_type == 'all':
    min_prevalence = 0.2
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


f_max_all = []


for species_name in species_list_to_plot:

    f_mean = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_mean_mapgd']
    f_mean = numpy.asarray(f_mean)
    predicted_f_mean = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_f_mean_mapgd']
    predicted_f_mean = numpy.asarray(predicted_f_mean)
    f_mean_error = numpy.absolute(f_mean - predicted_f_mean)/f_mean

    f_var = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_var_mapgd']
    f_var = numpy.asarray(f_var)
    predicted_f_var = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_f_var_mapgd']
    predicted_f_var = numpy.asarray(predicted_f_var)
    f_var_error = numpy.absolute(f_var - predicted_f_var)/f_var

    observed_prevalence = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
    observed_prevalence = numpy.asarray(observed_prevalence)
    predicted_prevalence = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd']
    predicted_prevalence = numpy.asarray(predicted_prevalence)
    prevalence_error = numpy.absolute(observed_prevalence - predicted_prevalence)/observed_prevalence

    f_max_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_max_mapgd']
    f_max_mapgd = numpy.asarray(f_max_mapgd)

    n_non_zero_f = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
    n_non_zero_f = numpy.asarray(n_non_zero_f)


    # plot error
    n_observatins = list(set(n_non_zero_f))
    n_observatins.sort()
    f_mean_mean_error_all = [numpy.mean(f_mean_error[n_non_zero_f==i]) for i in n_observatins]
    f_var_mean_error_all = [numpy.mean(f_var_error[n_non_zero_f==i]) for i in n_observatins]
    prevalence_mean_error_all = [numpy.mean(prevalence_error[n_non_zero_f==i]) for i in n_observatins]

    #ax_mean_error.plot(n_observatins, f_mean_mean_error_all, c=species_color_map[species_name], lw=1.5, alpha=0.8, linestyle='-', zorder=2)
    #ax_var_error.plot(n_observatins, f_var_mean_error_all, c=species_color_map[species_name], lw=1.5, alpha=0.8, linestyle='-', zorder=2)
    #ax_prevalence_error.plot(n_observatins, prevalence_mean_error_all, c=species_color_map[species_name], lw=1.5, alpha=0.8, linestyle='-', zorder=2)

    #ax_mean_error.scatter(f_max_mapgd[n_non_zero_f == 6], f_mean_error[n_non_zero_f == 6], alpha=0.2, s=10, c=species_color_map[species_name])
    ax_mean_error.scatter(f_max_mapgd[n_non_zero_f == 6], f_mean_error[n_non_zero_f == 6], alpha=0.2, s=10, c=species_color_map[species_name])

    #ax_var_error.scatter(f_max_mapgd[n_non_zero_f == 6], f_var_error[n_non_zero_f == 6], c=species_color_map[species_name], lw=1.5, alpha=0.2, linestyle='-', zorder=2)
    ax_prevalence_error.scatter(f_max_mapgd, prevalence_error, c=species_color_map[species_name], lw=1.5, alpha=0.8, linestyle='-', zorder=2)


    #ax_f_max_vs_prevalence_error.plot(f_max_mapgd, prevalence_error, c=species_color_map[species_name], lw=1.5, alpha=0.8, linestyle='-', zorder=2)
    # bin f_max and prevalence error
    f_max_log10 = numpy.log10(f_max_mapgd)
    hist, bin_edges = numpy.histogram(f_max_log10, density=True, bins=40)
    bins_mean = numpy.asarray([0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )])
    prevalence_mre_bins = []
    bins_mean_to_keep = []
    for i in range(0, len(bin_edges)-1):
        idx_i = (f_max_log10 > bin_edges[i]) & (f_max_log10 < bin_edges[i+1])
        if len(idx_i) > 5:
            prevalence_mre_bins.append(numpy.mean(prevalence_error[idx_i]))
            bins_mean_to_keep.append(bins_mean[i])
    bins_mean_to_keep = numpy.asarray(bins_mean_to_keep)
    prevalence_mre_bins = numpy.asarray(prevalence_mre_bins)


    ax_f_max_vs_prevalence_error.plot(10**bins_mean_to_keep, 10**prevalence_mre_bins, c=species_color_map[species_name], lw=1.5, alpha=0.8, linestyle='-', zorder=2)




    # filter for sites with more than one non zero observation
    f_mean = f_mean[n_non_zero_f > min_n_non_zero_f]
    predicted_f_mean = predicted_f_mean[n_non_zero_f > min_n_non_zero_f]
    f_var = f_var[n_non_zero_f > min_n_non_zero_f]
    predicted_f_var = predicted_f_var[n_non_zero_f > min_n_non_zero_f]


    # left column plots
    if len(predicted_f_mean) > n_sample:

        sample_idx = numpy.random.choice(len(predicted_f_mean), n_sample, replace=False)

        f_mean = f_mean[sample_idx]
        predicted_f_mean = predicted_f_mean[sample_idx]

        f_var = f_var[sample_idx]
        predicted_f_var = predicted_f_var[sample_idx]

        observed_prevalence = observed_prevalence[sample_idx]
        predicted_prevalence = predicted_prevalence[sample_idx]

        f_max_mapgd = f_max_mapgd[sample_idx]
        n_non_zero_f = n_non_zero_f[sample_idx]


    ax_mean.scatter(f_mean, predicted_f_mean, alpha=0.2, s=10, c=species_color_map[species_name])
    ax_var.scatter(f_var, predicted_f_var, alpha=0.2, s=10, c=species_color_map[species_name])
    ax_prevalence.scatter(observed_prevalence, predicted_prevalence, alpha=0.2, s=10, c=species_color_map[species_name])
    ax_f_max_vs_prevalence.scatter(f_max_mapgd, n_non_zero_f, alpha=0.2, s=10, c=species_color_map[species_name])





# relationship between error and f_max
#hist, bin_edges = numpy.histogram(f_max_log10, density=True, bins=40)

#bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]

#prevalence_predicted_mean = []
#relative_error_mean = []
#for i in range(0, len(bin_edges)-1):
#    idx_i = (f_max_log10 > bin_edges[i]) & (f_max_log10 < bin_edges[i+1])
#    prevalence_predicted_mean.append(numpy.mean(numpy.log10(predicted_prevalence[idx_i])))
#    relative_error_mean.append(numpy.mean(relative_error[idx_i]))



ax_mean.set_xscale('log', basex=10)
ax_mean.set_yscale('log', basey=10)
ax_mean.plot([0.00001,1],[0.00001,1], lw=3,ls='--',c='k',zorder=1)
ax_mean.set_xlim([0.0002,0.2])
ax_mean.set_ylim([0.0002,0.2])
ax_mean.set_xlabel('Observed mean SNV frequency', fontsize=11)
ax_mean.set_ylabel('Predicted mean SNV frequency', fontsize=11)


ax_var.set_xscale('log', basex=10)
ax_var.set_yscale('log', basey=10)
ax_var.plot([0.00001,1],[0.00001,1], lw=3,ls='--',c='k',zorder=1)
ax_var.set_xlim([0.00008,0.09])
ax_var.set_ylim([0.00008,0.09])
ax_var.set_xlabel('Observed vairance of SNV frequencies', fontsize=11)
ax_var.set_ylabel('Predicted vairance of SNV frequencies', fontsize=11)


ax_prevalence.set_xscale('log', basex=10)
ax_prevalence.set_yscale('log', basey=10)
ax_prevalence.plot([0.00001,1], [0.00001,1], lw=3,ls='--',c='k',zorder=1)
ax_prevalence.set_xlim([0.008,0.2])
ax_prevalence.set_ylim([0.0008,0.2])
ax_prevalence.set_xlabel('Observed SNV prevalence', fontsize=11)
ax_prevalence.set_ylabel('Predicted SNV prevalence', fontsize=11)


ax_f_max_vs_prevalence.set_xscale('log', basex=10)
ax_f_max_vs_prevalence.set_yscale('log', basey=10)
ax_f_max_vs_prevalence.set_xlim([0.008,1.04])
#ax_f_max_vs_prevalence.set_ylim([0.0008,0.2])
ax_f_max_vs_prevalence.set_xlabel('Max. SNV frequency within a host', fontsize=11)
ax_f_max_vs_prevalence.set_ylabel('Number of hosts where SNV is present', fontsize=11)




ax_mean_error.set_xscale('log', basex=10)
ax_mean_error.set_yscale('log', basey=10)
#ax_mean_error.set_ylim([0.09,1.1])
ax_mean_error.set_xlabel('Number of hosts where SNV is present', fontsize=11)
ax_mean_error.set_ylabel('MRE of predicted mean', fontsize=11)


ax_var_error.set_xscale('log', basex=10)
ax_var_error.set_yscale('log', basey=10)
#ax_var_error.set_ylim([0.09,1.1])
ax_var_error.set_xlabel('Number of hosts where SNV is present', fontsize=11)
ax_var_error.set_ylabel('MRE of predicted variance', fontsize=11)


ax_prevalence_error.set_xscale('log', basex=10)
ax_prevalence_error.set_yscale('log', basey=10)
ax_prevalence_error.set_xlabel('Number of hosts where SNV is present', fontsize=11)
ax_prevalence_error.set_ylabel('MRE of predicted prevalence', fontsize=11)


ax_f_max_vs_prevalence_error.set_xscale('log', basex=10)
ax_f_max_vs_prevalence_error.set_yscale('log', basey=10)
ax_f_max_vs_prevalence_error.set_xlabel('Max. SNV frequency within a host', fontsize=11)
ax_f_max_vs_prevalence_error.set_ylabel('MRE of predicted prevalence', fontsize=11)
#ax_f_max_vs_prevalence_error.set_ylim([0.001,110.1])




#ax_f_max_vs_prevalence.set_xlim([0.008,1.05])


fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.22)
fig.savefig("%sprediction_summary_%s_%s.png" % (config.analysis_directory, clade_type, variant_type), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
