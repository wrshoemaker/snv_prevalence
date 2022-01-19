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
import plot_utils
import stats_utils

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from scipy.stats import gamma, gaussian_kde, ks_2samp

import calculate_predicted_prevalence_mapgd



data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])


#clade_type = 'no_strains'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'
#max_n_occurances = 7

species_color_map, ordered_species_list = plot_utils.get_species_color_map()

prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()


species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()

n_species_row=5

nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]



gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
fig = plt.figure(figsize = (40, 20))

species_count = 0
for column_idx, column in enumerate(nested_species_list):

    for row_idx, row in enumerate(column):

        species_name = row

        print(row)

        f_mapgd_all = prevalence_dict_mapgd[species_name]['all'][pi_type][variant_type]['f_no_zeros_mapgd']
        f_mapgd_all = numpy.log10(f_mapgd_all)
        f_mapgd_all = (f_mapgd_all - numpy.mean(f_mapgd_all)) / numpy.std(f_mapgd_all)

        f_mapgd_no_strains = prevalence_dict_mapgd[species_name]['no_strains'][pi_type][variant_type]['f_no_zeros_mapgd']
        f_mapgd_no_strains = numpy.log10(f_mapgd_no_strains)
        f_mapgd_no_strains = (f_mapgd_no_strains - numpy.mean(f_mapgd_no_strains)) / numpy.std(f_mapgd_no_strains)


        hist_all, bin_edges_all = numpy.histogram(f_mapgd_all, density=True, bins=20)
        bins_mean_all = [0.5 * (bin_edges_all[i] + bin_edges_all[i+1]) for i in range(0, len(bin_edges_all)-1 )]

        hist_no_strains, bin_edges_no_strains = numpy.histogram(f_mapgd_no_strains, density=True, bins=20)
        bins_mean_no_strains = [0.5 * (bin_edges_no_strains[i] + bin_edges_no_strains[i+1]) for i in range(0, len(bin_edges_no_strains)-1 )]


        #ax.plot(bins_mean_all, hist_all, markersize=10, c='k', marker='o')
        #ax.plot(bins_mean_no_strains, hist_no_strains, markersize=10, c='k', marker='o', mfc='none', mew=2)

        #ax.scatter(bins_mean_all, hist_all, alpha=1, s=25, c=species_color_map[species_name])
        #ax.scatter(bins_mean_no_strains, hist_no_strains, alpha=1, s=25, linewidth=2, facecolor='none', edgecolor=species_color_map[species_name])
        freqs_concat = numpy.concatenate([f_mapgd_all, f_mapgd_no_strains])
        freq_range = numpy.linspace(min(freqs_concat)-0.2, max(freqs_concat)+0.2, num=100)


        survival_array_all = [sum(f_mapgd_all>=i)/len(f_mapgd_all) for i in freq_range]
        survival_array_all = numpy.asarray(survival_array_all)

        survival_array_no_strains = [sum(f_mapgd_no_strains>=i)/len(f_mapgd_no_strains) for i in freq_range]
        survival_array_no_strains = numpy.asarray(survival_array_no_strains)

        ax = fig.add_subplot(gs[column_idx, row_idx])

        ax.set_title(figure_utils.get_pretty_species_name(species_name), fontsize=12, fontweight='bold', color='k' )
        ax.plot(freq_range, survival_array_all, ls='-', lw=3, c='k', alpha=0.8, zorder=2, label='All')
        ax.plot(freq_range, survival_array_no_strains, ls='--', lw=3, c='k', alpha=0.8, zorder=2, label='No strains')

        #d, p = stats_utils.permutation_ks_test(f_mapgd_all, f_mapgd_no_strains, n_iter=10000)
        #print(d, p)

        if (row_idx == 0):
            ax.set_ylabel('Cumulative density', fontsize=14)


        if column_idx == len(nested_species_list)-1:
            ax.set_xlabel('Rescaled ' + r'$\mathrm{log}_{10}$' + ' SNV frequencies', fontsize=14)


        if column_idx == len(nested_species_list)-2:

            if row_idx > len(nested_species_list[-1])-1:

                ax.set_xlabel('Rescaled ' + r'$\mathrm{log}_{10}$' + 'SNV frequencies', fontsize=14)

        if (column_idx == 0) and (row_idx == 0):
            ax.legend(loc='lower left', prop={'size': 11})

        #species_count += 1





fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.25)
fig.savefig("%sfreq_dist_all_vs_no_strains_%s.png" % (config.analysis_directory, variant_type), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
