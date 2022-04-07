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

import calculate_predicted_prevalence

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])


prevalence_dict = calculate_predicted_prevalence.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict.keys())
species_list.sort()

n_species_row=5

nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]

def make_plot(clade_type, pi_type, variant_type='4D'):

    gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
    fig = plt.figure(figsize = (40, 20))

    for column_idx, column in enumerate(nested_species_list):

        for row_idx, row in enumerate(column):


                ax = fig.add_subplot(gs[column_idx, row_idx])


                predicted_prevalence = prevalence_dict[row][clade_type][pi_type][variant_type]['predicted_prevalence_to_plot']
                predicted_prevalence = numpy.asarray(predicted_prevalence)

                observed_prevalence = prevalence_dict[row][clade_type][pi_type][variant_type]['observed_prevalence_to_plot']
                observed_prevalence = numpy.asarray(observed_prevalence)

                all_ = numpy.concatenate([predicted_prevalence,observed_prevalence])


                # Calculate the point density
                xy = numpy.vstack([observed_prevalence, predicted_prevalence])
                z = gaussian_kde(xy)(xy)

                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = observed_prevalence[idx], predicted_prevalence[idx], z[idx]

                max_ = max(all_)*1.1
                min_ = min(all_)*0.8

                ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
                ax.set_xlim([min_, max_])
                ax.set_ylim([min_, max_])


                ax.scatter(x, y, c=z, cmap="Blues", s=30, alpha=0.5, edgecolor='', zorder=1)

                ax.set_xscale('log', basex=10)
                ax.set_yscale('log', basey=10)

                ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )


                if (row_idx == 0):
                    ax.set_ylabel('Predicted SNV prevalence', fontsize=12)


                if column_idx == len(nested_species_list)-1:
                    ax.set_xlabel('Observed SNV prevalence', fontsize=12)


                if column_idx == len(nested_species_list)-2:

                    if row_idx > len(nested_species_list[-1])-1:

                        ax.set_xlabel('Observed SNV prevalence', fontsize=12)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%spredicted_observed_prevalence_%s_%s_%s.png" % (config.analysis_directory, clade_type, pi_type, variant_type), format='png', bbox_inches = "tight", pad_inches = 0.4)
    plt.close()



clade_types = ['all','largest_clade', 'no_strains']
pi_types = ['pi_exclude_boundary', 'pi_include_boundary']

for allowed_variant_type in allowed_variant_types:
    for clade_type in clade_types:
        for pi_type in pi_types:
            make_plot(clade_type, pi_type, variant_type=variant_type)
