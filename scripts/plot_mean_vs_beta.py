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


#prevalence_dict = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()
prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()
species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()
n_species_row=5
nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]


def make_plot(clade_type, cov_cutoff):

    #clade_type = 'all'
    pi_type = 'pi_include_boundary'
    variant_type = '4D'

    #cov_cutoff = False



    gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
    fig = plt.figure(figsize = (40, 20))

    for column_idx, column in enumerate(nested_species_list):

        for row_idx, row in enumerate(column):

            ax = fig.add_subplot(gs[column_idx, row_idx])


            observed_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['f_mean_slm']
            observed_prevalence = numpy.asarray(observed_prevalence)

            predicted_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['f_beta_slm']
            predicted_prevalence = numpy.asarray(predicted_prevalence)


            predicted_prevalence_no_zeros = predicted_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]
            observed_prevalence_no_zeros = observed_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]


            #f_max_no_zeros = f_max[(observed_prevalence>0) & (predicted_prevalence>0)]
            if len(observed_prevalence_no_zeros) == 0:
                continue

            all_ = numpy.concatenate([predicted_prevalence_no_zeros,observed_prevalence_no_zeros])

            #re = numpy.absolute(observed_prevalence_no_zeros - predicted_prevalence_no_zeros) / observed_prevalence_no_zeros
            #re_all = numpy.absolute(observed_prevalence_all_no_zeros - predicted_prevalence_all_no_zeros) / observed_prevalence_all_no_zeros

            #mre = numpy.mean(re)
            #mre_all = numpy.mean(re_all)

            xy = numpy.vstack([observed_prevalence_no_zeros, predicted_prevalence_no_zeros])
            z = gaussian_kde(xy)(xy)
            # Sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            x, y, z = observed_prevalence_no_zeros[idx], predicted_prevalence_no_zeros[idx], z[idx]
            ax.scatter(x, y, c=z, cmap="Blues", s=90, alpha=0.9, edgecolor='', zorder=1)


            #ax.scatter(observed_prevalence, predicted_prevalence, c='dodgerblue', s=90, alpha=0.9, edgecolor='', zorder=1)

            max_ = max(all_)*1.1
            min_ = min(all_)*0.8

            #print(min(f_max_no_zeros))

            ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
            ax.set_xlim([min_, max_])
            ax.set_ylim([min_, max_])


            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)

            ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )


            if (row_idx == 0):
                ax.set_ylabel('Estimated ' + r'$\beta$', fontsize=12)


            if column_idx == len(nested_species_list)-1:
                ax.set_xlabel('Estimated mean frequency, ' + r'$\bar{f}$', fontsize=12)


            if column_idx == len(nested_species_list)-2:

                if row_idx > len(nested_species_list[-1])-1:

                    ax.set_xlabel('Estimated mean frequency, ' + r'$\bar{f}$', fontsize=12)


    #mre_all = numpy.asarray(mre_all)
    #print(sum(mre_all < 0.5) / len(mre_all))

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%sf_mean_vs_beta_%s_%s_%s.png" % (config.analysis_directory, clade_type, pi_type, variant_type), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
    plt.close()







#for clade_type in ['all', 'no_strains']:

make_plot('all', True)
