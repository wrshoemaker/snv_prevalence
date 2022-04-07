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

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_prevalence_mapgd

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

#clade_types = ['all','largest_clade', 'no_strains']

clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'


#prevalence_dict = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()
prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()
species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()
n_species_row=5
nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]


gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
fig = plt.figure(figsize = (40, 20))

for column_idx, column in enumerate(nested_species_list):

    for row_idx, row in enumerate(column):

        ax = fig.add_subplot(gs[column_idx, row_idx])

        predicted_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm']
        predicted_prevalence = numpy.asarray(predicted_prevalence)

        observed_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
        observed_prevalence = numpy.asarray(observed_prevalence)

        predicted_prevalence_no_zeros = predicted_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]
        observed_prevalence_no_zeros = observed_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]

        #f_max_no_zeros = f_max[(observed_prevalence>0) & (predicted_prevalence>0)]
        if len(observed_prevalence_no_zeros) == 0:
            continue

        predicted_prevalence_best = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm_best']
        predicted_prevalence_best = numpy.asarray(predicted_prevalence_best)

        observed_prevalence_best = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm_best']
        observed_prevalence_best = numpy.asarray(observed_prevalence_best)

        predicted_prevalence_best_no_zeros = predicted_prevalence_best[(observed_prevalence_best>0) & (predicted_prevalence_best>0) ]
        observed_prevalence_best_no_zeros = observed_prevalence_best[(observed_prevalence_best>0) & (predicted_prevalence_best>0) ]


        error = numpy.absolute(predicted_prevalence_no_zeros - observed_prevalence_no_zeros) / observed_prevalence_no_zeros
        error_best = numpy.absolute(predicted_prevalence_best_no_zeros - observed_prevalence_best_no_zeros) / observed_prevalence_best_no_zeros


        # fraction of SNVs where estimator improves prediction

        predicted_prevalence_best_no_zeros_all = predicted_prevalence_best[(observed_prevalence_best>0) & (predicted_prevalence_best>0) & (observed_prevalence>0) & (predicted_prevalence>0) ]
        observed_prevalence_best_no_zeros_all = observed_prevalence_best[(observed_prevalence_best>0) & (predicted_prevalence_best>0) & (observed_prevalence>0) & (predicted_prevalence>0) ]
        predicted_prevalence_no_zeros_all = predicted_prevalence[(observed_prevalence_best>0) & (predicted_prevalence_best>0) & (observed_prevalence>0) & (predicted_prevalence>0) ]
        observed_prevalence_no_zeros_all = observed_prevalence[(observed_prevalence_best>0) & (predicted_prevalence_best>0) & (observed_prevalence>0) & (predicted_prevalence>0) ]

        error_all = numpy.absolute(observed_prevalence_no_zeros_all-predicted_prevalence_no_zeros_all)/observed_prevalence_no_zeros_all
        error_best_all = numpy.absolute(observed_prevalence_best_no_zeros_all-predicted_prevalence_best_no_zeros_all)/observed_prevalence_best_no_zeros_all

        print(row, sum(error_best_all < error_all) / len(error_all))


        error_log10 = numpy.log10(error)
        error_best_log10 = numpy.log10(error_best)

        error_log10_no_nan = error_log10[numpy.isfinite(error_log10)]
        error_best_log10_no_nan = error_best_log10[numpy.isfinite(error_best_log10)]

        error_no_nan = 10**error_log10_no_nan
        error_best_no_nan = 10**error_best_log10_no_nan

        min_x = min([min(error_log10), min(error_best_log10)])
        max_x = max([max(error_log10), max(error_best_log10)])

        x_range = numpy.logspace(min_x, max_x, num=100, endpoint=True, base=10.0)


        survival_error = [sum(error_no_nan >= i)/len(error_no_nan) for i in x_range]
        survival_error_best = [sum(error_best_no_nan >= i)/len(error_best_no_nan) for i in x_range]

        survival_error = numpy.asarray(survival_error)
        survival_error_best = numpy.asarray(survival_error_best)

        ax.plot(x_range, survival_error, ls='-', lw=3, c='k', zorder=2, label='Gamma estimator')
        ax.plot(x_range, survival_error_best, ls='--', lw=3, c='k', zorder=2, label='Truncated gamma estimator')



        #hist, bins = numpy.histogram(error, bins=12)
        #logbins = numpy.logspace(numpy.log10(bins[0]), numpy.log10(bins[-1]), len(bins))
        #ax.hist(error, bins=logbins, color='k', density=True, linewidth=2, linestyle='-', histtype='step')

        #error_best_log10 = numpy.log10(error_best)
        #hist_best, bin_edges_best = numpy.histogram(error_best_log10, density=True, bins=12)
        #bins_mean_best = [0.5 * (bin_edges_best[i] + bin_edges_best[i+1]) for i in range(0, len(bin_edges_best)-1 )]
        #bins_mean_best = numpy.asarray(bins_mean_best)
        #hist_to_plot_best = hist_best[hist_best>0]
        #bins_mean_to_plot_best = bins_mean_best[hist_best>0]
        #ax.plot(10**bins_mean_to_plot_best, hist_to_plot_best, ls='--', c='k', lw=2)

        #hist_best, bins_best = numpy.histogram(error_best, bins=12)
        #logbins_best = numpy.logspace(numpy.log10(bins_best[0]), numpy.log10(bins_best[-1]),len(bins_best))
        #ax.hist(error_best, bins=logbins_best, color='k', density=True, linewidth=2, linestyle='--', histtype='step')



        #histtype='step', color='k', lw=3, alpha=0.8, bins= 20, density=True, zorder=2


        #ax.scatter(observed_prevalence, predicted_prevalence, c='dodgerblue', s=90, alpha=0.9, edgecolor='', zorder=1)

        #max_ = max(all_)*1.1
        #min_ = min(all_)*0.8

        #print(min(f_max_no_zeros))

        #ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
        #ax.set_xlim([min_, max_])
        #ax.set_ylim([min_, max_])


        ax.set_xscale('log', basex=10)
        #ax.set_yscale('log', basey=10)

        ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )


        if (row_idx == 0) and (column_idx==0):
            ax.legend(loc="lower left", fontsize=8)


        if (row_idx == 0):
            ax.set_ylabel('Fraction of SNVs ' + r'$\geq \epsilon$', fontsize=12)


        if column_idx == len(nested_species_list)-1:
            ax.set_xlabel('Relative error of SLM, ' + r'$\epsilon$', fontsize=12)


        if column_idx == len(nested_species_list)-2:

            if row_idx > len(nested_species_list[-1])-1:

                ax.set_xlabel('Relative error of SLM, ' + r'$\epsilon$', fontsize=12)


#mre_all = numpy.asarray(mre_all)
#print(sum(mre_all < 0.5) / len(mre_all))

fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# dpi = 600
fig.savefig("%spredicted_observed_prevalence_error_mapgd_slm_%s_%s.png" % (config.analysis_directory, clade_type, variant_type), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
plt.close()
