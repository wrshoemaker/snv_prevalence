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

        min_distance_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['min_distance_slm']
        min_distance_slm = numpy.asarray(min_distance_slm)

        # fraction of SNVs where estimator improves prediction
        idx_ = (observed_prevalence_best>0) & (predicted_prevalence_best>0) & (observed_prevalence>0) & (predicted_prevalence>0)

        predicted_prevalence_best_no_zeros_all = predicted_prevalence_best[idx_]
        observed_prevalence_best_no_zeros_all = observed_prevalence_best[idx_]
        predicted_prevalence_no_zeros_all = predicted_prevalence[idx_]
        observed_prevalence_no_zeros_all = observed_prevalence[idx_]

        min_distance_slm_no_zeros_all = min_distance_slm[idx_]

        error_all = numpy.absolute(observed_prevalence_no_zeros_all-predicted_prevalence_no_zeros_all)/observed_prevalence_no_zeros_all
        error_best_all = numpy.absolute(observed_prevalence_best_no_zeros_all-predicted_prevalence_best_no_zeros_all)/observed_prevalence_best_no_zeros_all

        #error_prevalence_improvement = (error_best_all-error_all)/error_all
        error_ratio = error_best_all/error_all


        min_distance_slm_no_zeros_all_better = min_distance_slm_no_zeros_all[error_ratio < 1]
        error_ratio_better = error_ratio[error_ratio < 1]

        min_distance_slm_no_zeros_all_worse = min_distance_slm_no_zeros_all[error_ratio >= 1]
        error_ratio_worse = error_ratio[error_ratio >= 1]


        all_ = numpy.concatenate([min_distance_slm_no_zeros_all, error_ratio])

        xy_better = numpy.vstack([min_distance_slm_no_zeros_all_better, error_ratio_better])
        z_better = gaussian_kde(xy_better)(xy_better)
        # Sort the points by density, so that the densest points are plotted last
        idx_better = z_better.argsort()
        x, y, z = min_distance_slm_no_zeros_all_better[idx_better], error_ratio_better[idx_better], z_better[idx_better]
        ax.scatter(min_distance_slm_no_zeros_all_better, error_ratio_better, c=z, cmap="Blues", s=30, alpha=0.9, edgecolor='', zorder=1)

        xy_worse = numpy.vstack([min_distance_slm_no_zeros_all_worse, error_ratio_worse])
        z_worse = gaussian_kde(xy_worse)(xy_worse)
        # Sort the points by density, so that the densest points are plotted last
        idx_worse = z_worse.argsort()
        x, y, z = min_distance_slm_no_zeros_all_worse[idx_worse], error_ratio_worse[idx_worse], z_worse[idx_worse]
        ax.scatter(min_distance_slm_no_zeros_all_worse, error_ratio_worse, c=z, cmap="Reds", s=30, alpha=0.9, edgecolor='', zorder=1)



        #max_ = max(all_)*1.1
        #min_ = min(all_)*0.8

        ax.axhline(1, lw=4, ls='--',color='k', zorder=2)

        #ax.set_xlim([min(min_distance_slm_no_zeros_all)*0.8, max(min_distance_slm_no_zeros_all)*1.1])
        ax.set_ylim([min(error_ratio)*0.8, max(error_ratio)*1.1])

        #ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)

        ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )


        #if (row_idx == 0) and (column_idx==0):
        #    ax.legend(loc="lower left", fontsize=8)


        if (row_idx == 0):
            ax.set_ylabel('Ratio of SLM error, ' + r'$\epsilon_{\mathrm{Trunc}} /  \epsilon_{\mathrm{Naive}}$', fontsize=12)


        if column_idx == len(nested_species_list)-1:
            ax.set_xlabel('Truncated estimator distance', fontsize=12)


        if column_idx == len(nested_species_list)-2:

            if row_idx > len(nested_species_list[-1])-1:

                ax.set_xlabel('Truncated estimator distance', fontsize=12)


#mre_all = numpy.asarray(mre_all)
#print(sum(mre_all < 0.5) / len(mre_all))

fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# dpi = 600
fig.savefig("%struncated_gamma_change_vs_distance_%s_%s.png" % (config.analysis_directory, clade_type, variant_type), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
plt.close()
