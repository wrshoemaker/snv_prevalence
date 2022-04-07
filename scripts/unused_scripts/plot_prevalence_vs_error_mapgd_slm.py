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

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_prevalence_mapgd

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()

species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()
n_species_row=5
nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]

clade_type = 'all'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
max_n_occurances = 7


gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
fig = plt.figure(figsize = (40, 20))

#mre_all = []

for column_idx, column in enumerate(nested_species_list):

    for row_idx, row in enumerate(column):

        ax = fig.add_subplot(gs[column_idx, row_idx])

        for variant_type in allowed_variant_types:

            predicted_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm']
            predicted_prevalence = numpy.asarray(predicted_prevalence)

            observed_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
            observed_prevalence = numpy.asarray(observed_prevalence)

            predicted_prevalence_no_zeros = predicted_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]
            observed_prevalence_no_zeros = observed_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]

            #f_max_no_zeros = f_max[(observed_prevalence>0) & (predicted_prevalence>0)]
            if len(observed_prevalence_no_zeros) == 0:
                continue

            all_ = numpy.concatenate([predicted_prevalence_no_zeros,observed_prevalence_no_zeros])

            #re = numpy.absolute(observed_prevalence_no_zeros - predicted_prevalence_no_zeros) / observed_prevalence_no_zeros
            #re_all = numpy.absolute(observed_prevalence_all_no_zeros - predicted_prevalence_all_no_zeros) / observed_prevalence_all_no_zeros

            relative_error = numpy.absolute(observed_prevalence_no_zeros - predicted_prevalence_no_zeros) / observed_prevalence_no_zeros

            observed_prevalence_no_zeros_log10 = numpy.log10(observed_prevalence_no_zeros)
            hist_all, bin_edges_all = numpy.histogram(observed_prevalence_no_zeros_log10, density=True, bins=20)
            bins_mean_all = [0.5 * (bin_edges_all[i] + bin_edges_all[i+1]) for i in range(0, len(bin_edges_all)-1 )]
            bins_mean_all_to_keep = []
            relative_error_mean_bins = []
            for i in range(0, len(bin_edges_all)-1 ):
                relative_error_i = relative_error[ (observed_prevalence_no_zeros_log10>=bin_edges_all[i]) & (observed_prevalence_no_zeros_log10<bin_edges_all[i+1])]
                if len(relative_error_i) < 10:
                    continue
                bins_mean_all_to_keep.append(bins_mean_all[i])
                relative_error_mean_bins.append(numpy.mean(relative_error_i))

            bins_mean_all_to_keep = numpy.asarray(bins_mean_all_to_keep)
            relative_error_mean_bins = numpy.asarray(relative_error_mean_bins)
            #xy = numpy.vstack([observed_prevalence_no_zeros, relative_error])
            #z = gaussian_kde(xy)(xy)
            # Sort the points by density, so that the densest points are plotted last
            #idx = z.argsort()
            #x, y, z = observed_prevalence_no_zeros[idx], relative_error[idx], z[idx]
            #ax.scatter(x, y, c=z, cmap="Blues", s=90, alpha=0.9, edgecolor='', zorder=1)

            ax.scatter(10**bins_mean_all_to_keep, relative_error_mean_bins, c=prevalence_utils.variant_color_dict[variant_type], s=90, alpha=0.9, edgecolor='', zorder=1)

            max_ = max(all_)*1.1
            min_ = min(all_)*0.8

            #ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
            ax.set_xlim([min_, max_])
            ax.set_ylim([0, max(relative_error_mean_bins)])
            ax.set_xscale('log', basex=10)
            #ax.set_yscale('log', basey=10)
            ax.axhline(y=0, color='k', linestyle=':', lw = 3, zorder=1)
            ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=14, fontweight='bold', color='k')


            if (row_idx == 0):
                ax.set_ylabel('Relative error', fontsize=14)

            if column_idx == len(nested_species_list)-1:
                ax.set_xlabel('Observed allele prevalence', fontsize=14)

            if column_idx == len(nested_species_list)-2:
                if row_idx > len(nested_species_list[-1])-1:
                    ax.set_xlabel('Observed allele prevalence', fontsize=14)


fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# dpi = 600
fig.savefig("%sprevalence_vs_error_mapgd_slm.pdf" % (config.analysis_directory), format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
