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

import calculate_predicted_prevalence_mapgd
import calculate_predicted_prevalence

import plot_utils


#plt.rcParams["text.usetex"] =True


pi_type = 'pi_include_boundary'
variant_type = '4D'
#clade_type = 'no_strains'

clade_types = ['all', 'no_strains']

species_color_map, ordered_species_list = plot_utils.get_species_color_map()

prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()


ls_dict = {'all': '--', 'no_strains':':'}

max_n_occurances = 7
range_max_n_occurances = range(1, max_n_occurances+1)


#fig, ax = plt.subplots(figsize=(4,4))

fig = plt.figure(figsize = (4.5, 12)) #
fig.subplots_adjust(bottom= 0.15)

ax_mean = plt.subplot2grid((3, 1), (0,0))
ax_var = plt.subplot2grid((3, 1), (1,0))
ax_prevalence = plt.subplot2grid((3, 1), (2,0))

axes_all = [ax_mean, ax_var, ax_prevalence]

predictions = ['mean', 'var', 'prevalence']


error_dict = {}
for p in predictions:
    error_dict[p] = {}
    for clade_type in clade_types:
        error_dict[p][clade_type] = {}
        for i in range_max_n_occurances:
            error_dict[p][clade_type][i] = {}
            error_dict[p][clade_type][i] = []
            #error_dict[clade_type][i]['lower_quartile'] = []
            #error_dict[clade_type][i]['upper_quartile'] = []


error_difference_dict = {}
for p in predictions:
    error_difference_dict[p] = {}
    for species_name in species_list:
        error_difference_dict[p][species_name] = {}
        for i in range_max_n_occurances:
            error_difference_dict[p][species_name][i] = {}
            #for clade_type in clade_types:
            #    error_difference_dict[species_name][i]



species_list_no_strains = []
for species_name in species_list:

    #error_dict[species_name] = {}

    for clade_type in clade_types:

        #error_dict[species_name][clade_type] = {}

        n_non_zero_frequencies = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
        n_non_zero_frequencies = numpy.asarray(n_non_zero_frequencies)

        observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
        observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)
        predicted_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd']
        predicted_prevalence_mapgd = numpy.asarray(predicted_prevalence_mapgd)
        rel_error_mapgd = numpy.absolute(observed_prevalence_mapgd - predicted_prevalence_mapgd) / observed_prevalence_mapgd


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



        n_non_zero_frequencies_counts = list(set(n_non_zero_frequencies))
        n_non_zero_frequencies_counts.sort()

        #mean_rel_error = [numpy.median(rel_error_mapgd[n_non_zero_frequencies==i]) for i in  n_non_zero_frequencies_counts]
        #median_all = []
        #n_non_zero_frequencies_counts_all = []
        for i in n_non_zero_frequencies_counts:
            if i not in range_max_n_occurances:
                continue
            #if i not in error_dict[p][clade_type]:
            #    continue
            prevalence_rel_error_i = rel_error_mapgd[n_non_zero_frequencies==i]
            f_mean_rel_error_i = f_mean_error[n_non_zero_frequencies==i]
            f_var_rel_error_i = f_var_error[n_non_zero_frequencies==i]


            prevalence_median = numpy.median(prevalence_rel_error_i)
            prevalence_mean = numpy.mean(prevalence_rel_error_i)

            f_mean_median = numpy.median(f_mean_rel_error_i)
            f_mean_mean = numpy.mean(f_mean_rel_error_i)

            f_var_median = numpy.median(f_var_rel_error_i)
            f_var_mean = numpy.mean(f_var_rel_error_i)



            #lower_quartile = numpy.quantile(rel_error_i, q=0.25)
            #upper_quartile = numpy.quantile(rel_error_i, q=0.75)

            #error_dict['prevalence'][clade_type][i].append(prevalence_median)
            #error_dict['mean'][clade_type][i].append(f_mean_median)
            #error_dict['var'][clade_type][i].append(f_var_median)

            #n_non_zero_frequencies_counts_all.append(i)
            #median_all.append(median)

            error_difference_dict['prevalence'][species_name][i][clade_type] = prevalence_mean
            error_difference_dict['mean'][species_name][i][clade_type] = f_mean_mean
            error_difference_dict['var'][species_name][i][clade_type] = f_var_mean

            #error_dict[clade_type][i]['median'].append(median)
            #error_dict[clade_type][i]['lower_quartile'].append(lower_quartile)
            #error_dict[clade_type][i]['upper_quartile'].append(upper_quartile)


        #error_dict[species_name][clade_type]['n_non_zero_frequencies_counts'] = numpy.asarray(n_non_zero_frequencies_counts)
        #error_dict[species_name][clade_type]['mean_rel_error'] = mean_rel_error

        #ax.plot(n_non_zero_frequencies_counts_all, median_all, c='k', lw=2.5, alpha=0.6, linestyle=ls_dict[clade_type], zorder=2)




#for clade_type in clade_types:
#    median_rel_error_all = [numpy.median(error_dict[clade_type][i]) for i in range_max_n_occurances]
#    lower_quartile_rel_error_all = [numpy.quantile(error_dict[clade_type][i], q=0.25) for i in range_max_n_occurances]
#    upper_quartile_rel_error_all = [numpy.quantile(error_dict[clade_type][i], q=0.75) for i in range_max_n_occurances]

#    #rint(median_rel_error_all, lower_quartile_rel_error_all, upper_quartile_rel_error_all)
#    #ax.plot(range_max_n_occurances, median_rel_error_all, c='k', lw=2.5, alpha=0.6, linestyle=ls_dict[clade_type], zorder=2)

#    #ax.errorbar(range_max_n_occurances, median_rel_error_all, yerr = [lower_quartile_rel_error_all, upper_quartile_rel_error_all], fmt = 'o', alpha = 0.5, \
#    #    barsabove = True, marker = '.', mfc = 'k', mec = 'k', c = 'k', zorder=1)
#    #plt.scatter(time_points_set, Ls, c='#175ac6', marker = 'o', s = 70, \
#    #    edgecolors='#244162', linewidth = 0.6, alpha = 0.5, zorder=2)#, edgecolors='none')


for p_idx, p in enumerate(predictions):

    ax = axes_all[p_idx]

    for species_name in species_list:
        range_max_n_occurances_to_plot = [i for i in range_max_n_occurances if ('no_strains' in error_difference_dict[p][species_name][i])]
        to_plot = [error_difference_dict[p][species_name][i]['no_strains']/error_difference_dict[p][species_name][i]['all'] for i in range_max_n_occurances if ('no_strains' in error_difference_dict[p][species_name][i])]

        ax.plot(range_max_n_occurances_to_plot, to_plot, c=species_color_map[species_name], lw=1.5, alpha=0.9, linestyle='-', zorder=2)


    #ax.plot(n_non_zero_frequencies_counts, mean_rel_error, c='k', lw=2.5, alpha=0.6, linestyle=ls_dict[clade_type], zorder=2)

    ax.axhline(y=1, color='k', linestyle=':', lw = 3, label='Equal MRE', zorder=1)

    ax.set_xlim([1, max_n_occurances])
    ax.set_yscale('log', basey=10)

    if p == 'mean':
        y_legend = 'Change in mean MRE due to strains, ' + r'$\epsilon_{\mathrm{No \, strains}} / \epsilon_{\mathrm{All}}$'
    elif p == 'var':
        y_legend = 'Change in variance MRE due to strains, ' + r'$\epsilon_{\mathrm{No \, strains}} / \epsilon_{\mathrm{All}}$'
    else:
        y_legend = 'Change in prevalence MRE due to strains, ' + r'$\epsilon_{\mathrm{No \, strains}} / \epsilon_{\mathrm{All}}$'


    ax.set_xlabel('Number of hosts where a SNV is present', fontsize=10)
    ax.set_ylabel(y_legend, fontsize=10)

    ax.legend(loc="lower right", fontsize=8)


fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.22)
fig.savefig("%serror_comparison.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
