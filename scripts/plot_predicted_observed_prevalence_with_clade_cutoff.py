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

import calculate_predicted_occupancy

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])
clade_types = ['all','largest_clade']

clade_type_label_dict = {'all': 'All hosts', 'largest_clade': 'Largest clade'}


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict.keys())
species_list.sort()

n_species_row=4

nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]

gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
fig = plt.figure(figsize = (26, 20))

for column_idx, column in enumerate(nested_species_list):

    for row_idx, row in enumerate(column):

        for clade_type_idx, clade_type in enumerate(clade_types):

            ax_ = fig.add_subplot(gs[column_idx, (2*row_idx) + clade_type_idx])

            predicted_prevalence = prevalence_dict[row][clade_type]['predicted_prevalence_to_plot']
            predicted_prevalence = numpy.asarray(predicted_prevalence)

            observed_prevalence = prevalence_dict[row][clade_type]['observed_prevalence_to_plot']
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

            ax_.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
            ax_.set_xlim([min_, max_])
            ax_.set_ylim([min_, max_])


            ax_.scatter(x, y, c=z, cmap="Blues", s=30, alpha=0.5, edgecolor='', zorder=1)

            ax_.set_xscale('log', basex=10)
            ax_.set_yscale('log', basey=10)

            if clade_type_idx == 0:

                ax_.text(1.05, 1.06, figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k', ha='center', va='center', transform=ax_.transAxes  )


            if column_idx == 0:
                #ax_species.set_xlabel('Observed SNV prevalence', fontsize=12)
                ax_.text(0.5, 1.14, clade_type_label_dict[clade_type], fontsize=12, fontweight='bold', color='k', ha='center', va='center', transform=ax_.transAxes  )



            if (row_idx == 0) and (clade_type == 'all'):
                ax_.set_ylabel('Predicted SNV prevalence', fontsize=12)


            if column_idx == len(nested_species_list)-1:
                ax_.set_xlabel('Observed SNV prevalence', fontsize=12)


            if column_idx == len(nested_species_list)-2:

                if row_idx > len(nested_species_list[-1])-1:

                    ax_.set_xlabel('Observed SNV prevalence', fontsize=12)








        #all_axis = fig.add_subplot(gs[column_idx, 2*row_idx])
        #largest_clade_axis = fig.add_subplot(gs[column_idx, (2*row_idx)+1])



        #predicted_prevalence = prevalence_dict[row]['predicted_observed_prevalence']['predicted_prevalence_to_plot']
        #predicted_prevalence = numpy.asarray(predicted_prevalence)

        #observed_prevalence = predicted_prevalence_dict[species]['predicted_observed_prevalence']['observed_prevalence_to_plot']
        #observed_prevalence = numpy.asarray(observed_prevalence)

        # y axis labels
        #if column_idx == 0:
        #    #ax_species.set_xlabel('Observed SNV prevalence', fontsize=12)
        #    all_axis.text(0.5, 1.15, "All hosts", fontsize=12, fontweight='bold', color='k', ha='center', va='center', transform=all_axis.transAxes  )
        #    largest_clade_axis.text(0.5, 1.15, "Largest clade", fontsize=12, fontweight='bold', color='k', ha='center', va='center', transform=largest_clade_axis.transAxes  )


        #if row_idx == 0:
        #    all_axis.set_ylabel('Predicted SNV prevalence', fontsize=12)


        #if column_idx == len(nested_species_list)-1:

        #    all_axis.set_xlabel('Observed SNV prevalence', fontsize=12)
        #    largest_clade_axis.set_xlabel('Observed SNV prevalence', fontsize=12)


        #if column_idx == len(nested_species_list)-2:

        #    if row_idx > len(nested_species_list[-1])-1:

        #        all_axis.set_xlabel('Observed SNV prevalence', fontsize=12)
        #        largest_clade_axis.set_xlabel('Observed SNV prevalence', fontsize=12)












fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
fig.savefig("%spredicted_observed_prevalence.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
