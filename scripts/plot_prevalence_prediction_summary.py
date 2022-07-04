from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path
from collections import Counter

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

numpy.random.seed(123456789)

#species_color_map, ordered_species_list = plot_utils.get_species_color_map()
species_color_map = plot_utils.species_color_map_genus

prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()#test=True)

species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()

clade_type = 'all'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
#variant_type = '4D'
best = True


max_n_occurances = 20
max_n_occurances_range = list(range(1, max_n_occurances+1))
iter=10000

def calculate_null_prevalence_vs_error():

    mean_all_dict = {}
    for j in max_n_occurances_range:
        mean_all_dict[j] = {}
        for species_name in species_list:
            mean_all_dict[j][species_name] = []

    for species_name in species_list:

        observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
        observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)

        predicted_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm']
        predicted_prevalence_mapgd = numpy.asarray(predicted_prevalence_mapgd)

        n_non_zero_f = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
        n_non_zero_f = numpy.asarray(n_non_zero_f)


        idx_ = (observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0) & (n_non_zero_f<=max_n_occurances)
        observed_prevalence_mapgd = observed_prevalence_mapgd[idx_]
        predicted_prevalence_mapgd = predicted_prevalence_mapgd[idx_]
        n_non_zero_f = n_non_zero_f[idx_]

        error = numpy.absolute(observed_prevalence_mapgd-predicted_prevalence_mapgd)/observed_prevalence_mapgd

        n_non_zero_f_dict = dict(Counter(n_non_zero_f.tolist()))
        n_non_zero_f_dict_items = n_non_zero_f_dict.items()
        n_non_zero_f_sorted = sorted(n_non_zero_f_dict_items)
        n_non_zero_f_sorted_n = [s[0] for s in n_non_zero_f_sorted]
        n_non_zero_f_sorted_perm_blocks = numpy.asarray([s[1] for s in n_non_zero_f_sorted])
        permutation_idx = numpy.asarray([sum(n_non_zero_f_sorted_perm_blocks[:block_idx+1]) for block_idx in range(len(n_non_zero_f_sorted)) ])

        for i in range(iter):
            error_permute = numpy.random.permutation(error)
            mean_error_blocks = [numpy.mean(s) for s in numpy.split(error_permute, permutation_idx)[:-1]]

            for n_non_zero_f_sorted_n_i_idx, n_non_zero_f_sorted_n_i in enumerate(n_non_zero_f_sorted_n):

                mean_all_dict[n_non_zero_f_sorted_n_i][species_name].append(mean_error_blocks[n_non_zero_f_sorted_n_i_idx])


    species_ci_dict = {}
    species_pooled_ci_dict = {}
    # calculate CIs for each species
    for n_non_zero_f in mean_all_dict.keys():

        error_species_all = []

        for species_name in mean_all_dict[n_non_zero_f].keys():

            error_species = mean_all_dict[n_non_zero_f][species_name]

            if len(error_species) == 0:
                continue

            error_species_all.extend(error_species)

            error_species = numpy.asarray(error_species)
            error_species = numpy.sort(error_species)

            #print(n_non_zero_f,species_name, len(error_species))

            error_species_ci_025 = error_species[int(len(error_species)*0.025)]
            error_species_ci_975 = error_species[int(len(error_species)*0.975)]

            if species_name not in species_ci_dict:
                species_ci_dict[species_name] = {}

            species_ci_dict[species_name][n_non_zero_f] = {}

            species_ci_dict[species_name][n_non_zero_f]['ci_025'] = error_species_ci_025
            species_ci_dict[species_name][n_non_zero_f]['ci_975'] = error_species_ci_975

        error_species_all = numpy.asarray(error_species_all)
        error_species_all = numpy.sort(error_species_all)
        error_species_all_ci_025 = error_species_all[int(len(error_species_all)*0.025)]
        error_species_all_ci_975 = error_species_all[int(len(error_species_all)*0.975)]

        species_pooled_ci_dict[n_non_zero_f] = {}
        species_pooled_ci_dict[n_non_zero_f]['ci_025'] = error_species_all_ci_025
        species_pooled_ci_dict[n_non_zero_f]['ci_975'] = error_species_all_ci_975

    return species_pooled_ci_dict, species_ci_dict



def make_plot(variant_type):

    fig = plt.figure(figsize = (8, 8))


    ax_error_dist = plt.subplot2grid((6, 6), (0, 0), colspan=3, rowspan=3)
    ax_obs_vs_pred = plt.subplot2grid((6, 6), (0, 3), colspan=3, rowspan=3)
    #ax_prevalence_vs_error = plt.subplot2grid((2, 2), (1, 0), colspan=1)
    ax_slope = plt.subplot2grid((6, 6), (3, 1), colspan=5, rowspan=3)

    #ax_all = [ax_obs_vs_pred, ax_error_dist, ax_slope]
    ax_error_dist.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_error_dist.transAxes)
    ax_obs_vs_pred.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_obs_vs_pred.transAxes)
    ax_slope.text(-0.29, 1, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_slope.transAxes)


    all_values = []
    slope_dict = {}
    for species_name in species_list:

        observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
        observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)

        predicted_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm']
        predicted_prevalence_mapgd = numpy.asarray(predicted_prevalence_mapgd)

        n_non_zero_f = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
        n_non_zero_f = numpy.asarray(n_non_zero_f)

        f_max_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_max_mapgd']
        f_max_mapgd = numpy.asarray(f_max_mapgd)



        idx_ = (observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0) & (f_max_mapgd<1)
        #idx_ = (observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0)
        observed_prevalence_mapgd = observed_prevalence_mapgd[idx_]
        predicted_prevalence_mapgd = predicted_prevalence_mapgd[idx_]
        n_non_zero_f = n_non_zero_f[idx_]

        error = numpy.absolute(observed_prevalence_mapgd-predicted_prevalence_mapgd)/observed_prevalence_mapgd

        print(species_name, sum(error<0.1)/len(error))

        observed_prevalence_mapgd_log10 = numpy.log10(observed_prevalence_mapgd)
        predicted_prevalence_mapgd_log10 = numpy.log10(predicted_prevalence_mapgd)

        # relationship between f_mean and prevalence
        hist, bin_edges = numpy.histogram(observed_prevalence_mapgd_log10, density=True, bins=10)
        bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]
        predicted_prevalence_mapgd_log10_mean = []
        mean_error = []
        bins_mean_to_plot = []
        for i in range(0, len(bin_edges)-1):
            idx_i = (observed_prevalence_mapgd_log10 > bin_edges[i]) & (observed_prevalence_mapgd_log10 <= bin_edges[i+1])
            if len(idx_i) <10:
                continue
            predicted_prevalence_mapgd_log10_mean.append(numpy.mean(predicted_prevalence_mapgd_log10[idx_i]))
            mean_error.append(numpy.mean(error[idx_i]))
            bins_mean_to_plot.append(bin_edges[i])

        all_values.extend(bins_mean_to_plot)
        all_values.extend(predicted_prevalence_mapgd_log10_mean)

        bins_mean_to_plot = numpy.asarray(bins_mean_to_plot)
        predicted_prevalence_mapgd_log10_mean = numpy.asarray(predicted_prevalence_mapgd_log10_mean)
        mean_error = numpy.asarray(mean_error)
        ax_obs_vs_pred.scatter(10**bins_mean_to_plot, 10**predicted_prevalence_mapgd_log10_mean, alpha=0.9, s=10, c=species_color_map[species_name], zorder=2)
        # label=figure_utils.get_pretty_species_name(species_name)


        error_log_rescaled = (numpy.log10(error) - numpy.mean(numpy.log10(error)))/numpy.std(numpy.log10(error))
        hist, bin_edges = numpy.histogram(error_log_rescaled, density=True, bins=20)
        bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]
        ax_error_dist.scatter(bins_mean, hist, alpha=1, s=10, c=species_color_map[species_name])


        #ax_error_dist.plot(10**bins_mean_to_plot, mean_error, lw=2, ls='-', c=species_color_map[species_name], alpha=0.5)
        #ax_error_dist.scatter(10**bins_mean_to_plot, mean_error, alpha=1, s=10, c=species_color_map[species_name])

        # plot observed prevalence vs. mean error
        #n_non_zero_f_set = list(set(n_non_zero_f.tolist()) & set(max_n_occurances_range))
        #mean_error_slm_n = [numpy.mean(error[n_non_zero_f == n]) for n in n_non_zero_f_set]
        #ax_prevalence_vs_error.plot(n_non_zero_f_set, mean_error_slm_n, lw=2, ls='-', marker='o', c=species_color_map[species_name], alpha=0.9)
        #ax_prevalence_vs_error.scatter(10**bins_mean_to_plot, mean_error, alpha=0.9, s=10, c=species_color_map[species_name], zorder=2)
        #ax_prevalence_vs_error.plot(10**bins_mean_to_plot, mean_error, lw=2, ls='-', marker='o', c=species_color_map[species_name], alpha=0.9)

        # relationship between prevalence and error
        #rho = numpy.corrcoef(observed_prevalence_mapgd_log10, error)[0,1]
        #slope, intercept, r_value, p_value, std_err = stats.linregress(observed_prevalence_mapgd_log10, error)
        error_log10 = numpy.log10(error)
        rho = numpy.corrcoef(observed_prevalence_mapgd_log10, error_log10)[0,1]
        #x_log10_range =  numpy.linspace(min(observed_prevalence_mapgd_log10) , max(observed_prevalence_mapgd_log10) , 10000)
        #y_fit_range = slope*x_log10_range + intercept
        #ax_prevalence_vs_error.plot(10**x_log10_range, y_fit_range, c=species_color_map[species_name], lw=2.5, linestyle='--', zorder=2)

        slope_null_all = []
        for i in range(iter):
            error_permute = numpy.random.permutation(error_log10)
            rho_null = numpy.corrcoef(observed_prevalence_mapgd_log10, error_permute)[0,1]
            #rho_null = numpy.corrcoef(observed_prevalence_mapgd_log10, error_permute)[0,1]
            #slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(observed_prevalence_mapgd_log10, error_permute)
            slope_null_all.append(rho_null)
        slope_null_all = numpy.asarray(slope_null_all)
        slope_null_all = numpy.sort(slope_null_all)
        median = slope_null_all[int(iter*0.5)]
        lower_ci = slope_null_all[int(iter*0.025)]
        upper_ci = slope_null_all[int(iter*0.975)]

        slope_dict[species_name] = {}
        slope_dict[species_name]['slope'] = rho
        slope_dict[species_name]['median'] = median
        slope_dict[species_name]['lower_ci'] = lower_ci
        slope_dict[species_name]['upper_ci'] = upper_ci


    min_ = 10**min(all_values)
    max_ = 1.06
    ax_obs_vs_pred.set_xlim([min_, max_])
    ax_obs_vs_pred.set_ylim([min_, max_])
    ax_obs_vs_pred.set_xscale('log', basex=10)
    ax_obs_vs_pred.set_yscale('log', basey=10)
    ax_obs_vs_pred.plot([min_, max_],[min_, max_], lw=3,ls='--',c='k',zorder=1, label='1:1')
    ax_obs_vs_pred.legend(loc='upper left', fontsize=11)
    ax_obs_vs_pred.set_xlabel('Observed prevalence', fontsize=11)
    ax_obs_vs_pred.set_ylabel('Predicted prevalence, SLM', fontsize=11)

    ax_obs_vs_pred.xaxis.set_tick_params(labelsize=8)
    ax_obs_vs_pred.yaxis.set_tick_params(labelsize=8)

    ax_error_dist.set_xlabel('Rescaled ' + r'$\mathrm{log}_{10}$' +  ' relative error', fontsize=11)
    ax_error_dist.set_ylabel('Probability density', fontsize=11)
    #ax_error_dist.set_xlabel('Observed prevalence', fontsize=11)
    #ax_error_dist.set_ylabel('Mean relative error', fontsize=11)


    #
    #ax_f_mean.set_ylim([-0.02, 0.9])
    #ax_error_dist.set_xlim([0.005, 1])
    #ax_error_dist.set_ylim([0.08, 1.05])
    #ax_error_dist.set_xscale('log', basex=10)
    ax_error_dist.set_xlim([-9, 3])
    ax_error_dist.set_ylim([0.00003, 4])
    ax_error_dist.set_yscale('log', basey=10)
    ax_error_dist.xaxis.set_tick_params(labelsize=8)
    ax_error_dist.yaxis.set_tick_params(labelsize=8)




    #species_pooled_ci_dict, species_ci_dict = calculate_null_prevalence_vs_error()

    #n_non_zero_f = list(species_pooled_ci_dict.keys())
    #n_non_zero_f.sort()
    #lower_ci = [species_pooled_ci_dict[n]['ci_025'] for n in n_non_zero_f]
    #upper_ci = [species_pooled_ci_dict[n]['ci_975'] for n in n_non_zero_f]

    #ax_prevalence_vs_error.fill_between(n_non_zero_f, lower_ci, upper_ci, alpha=0.3, color='grey', zorder=1, label='95% CI')
    #ax_prevalence_vs_error.plot(n_non_zero_f, lower_ci, color ='k', lw=1.5, zorder=1)
    #ax_prevalence_vs_error.plot(n_non_zero_f, upper_ci, color ='k', lw=1.5, zorder=1 )
    #ax_prevalence_vs_error.legend(loc='upper right', fontsize=11)

    #ax_prevalence_vs_error.set_xlabel('Observed prevalence', fontsize=12)
    #ax_prevalence_vs_error.set_ylabel('Relative error, SLM', fontsize=12)
    #ax_prevalence_vs_error.set_xscale('log', basex=10)
    #ax_prevalence_vs_error.set_xlim([0.005, max_])
    #ax_prevalence_vs_error.set_ylim([0, 0.9])




    # plot slope and null
    slope_tuple_list = [(species_name, slope_dict[species_name]['slope']) for species_name in slope_dict.keys()]
    slope_tuple_list.sort(key=lambda tup: tup[1])
    species_sorted = [s[0] for s in slope_tuple_list]
    species_sorted = species_sorted[::-1]
    species_sorted_pretty = [figure_utils.get_pretty_species_name(species_name) for species_name in species_sorted]

    for species_name_idx, species_name in enumerate(species_sorted):

        slope = slope_dict[species_name]['slope']
        median = slope_dict[species_name]['median']
        lower_ci = slope_dict[species_name]['lower_ci']
        upper_ci = slope_dict[species_name]['upper_ci']

        ax_slope.errorbar(median, species_name_idx, xerr=numpy.absolute(upper_ci), linestyle='-', markersize=5.5, marker='o', c='k', linewidth=2.5, alpha=1, zorder=1)
        ax_slope.scatter(slope, species_name_idx, alpha=1, s=55, c=species_color_map[species_name], zorder=2)

        print(slope)

    ax_slope.set_xlim([-0.8, 0.3])
    ax_slope.set_ylim([-0.5, len(species_sorted)])
    #ax_slope.axvline(x=0, color='k', linestyle='--', lw = 2, zorder=1)

    ax_slope.errorbar(-50000, -50000, xerr=10, linestyle='-', markersize=5, marker='o', c='k', linewidth=2.5, alpha=1, zorder=1, label='Null 95% CI')
    ax_slope.legend(loc='lower left', fontsize=11)

    ax_slope.set_xlabel('Correlation between observed ' +  r'$\mathrm{log}_{10}$' + ' prevalence and relative error', fontsize=11.5 )
    ax_slope.set_yticks(range(len(species_sorted)))
    ax_slope.set_yticklabels(species_sorted_pretty, fontsize=9)#, rotation = 45)


    #ax_slope.xaxis.set_tick_params(labelsize=8)
    ax_slope.xaxis.set_tick_params(labelsize=9)

    #ax_slopes.errorbar(list(range(1, 13)), slopes, slops_CIs,linestyle='-', marker='o', c='k', elinewidth=1.5, alpha=1, zorder=2)


    #fig.tight_layout()
    fig.subplots_adjust(wspace=1.25, hspace=0.8)
    fig.savefig("%sprevalence_summary_%s.pdf" % (config.analysis_directory, variant_type), format='pdf', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()





if __name__=='__main__':

    for variant_type in ['4D', '1D']:

        make_plot(variant_type)
