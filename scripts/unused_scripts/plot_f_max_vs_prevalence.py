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

n_species_row=5

nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]

gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
fig = plt.figure(figsize = (40, 20))

for column_idx, column in enumerate(nested_species_list):

    for row_idx, row in enumerate(column):

        ax = fig.add_subplot(gs[column_idx, row_idx])

        f_max = prevalence_dict[row]['all']['4D']['f_max']['to_plot']
        f_max = numpy.asarray(f_max)

        observed_prevalence = prevalence_dict[row]['all']['4D']['observed_prevalence_to_plot']
        observed_prevalence = numpy.asarray(observed_prevalence)

        all_ = numpy.concatenate([f_max,observed_prevalence])


        # Calculate the point density
        xy = numpy.vstack([f_max, observed_prevalence])
        z = gaussian_kde(xy)(xy)

        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = f_max[idx], observed_prevalence[idx], z[idx]

        #max_ = max(all_)*1.1
        #min_ = min(all_)*0.8

        #ax.set_xlim([min(f_max)*0.8, max(f_max)*1.1])
        ax.set_xlim([min(f_max)*0.8, 0.2])

        ax.set_ylim([min(observed_prevalence)*0.8, max(observed_prevalence)*1.1])


        ax.scatter(x, y, c=z, cmap="Blues", s=25, alpha=0.9, edgecolor='', zorder=1)
        #ax.scatter(f_max, observed_prevalence, s=25, alpha=0.9, edgecolor='', zorder=1)

        f_max_line = prevalence_dict[row]['all']['4D']['f_max_line']
        predicted_prevalence_line = prevalence_dict[row]['all']['4D']['f_max_vs_predicted_prevalence_line']

        f_max_line = numpy.asarray(f_max_line)
        predicted_prevalence_line = numpy.asarray(predicted_prevalence_line)

        ax.plot(10**f_max_line, 10**predicted_prevalence_line, c='k', lw = 5, ls ='--', label = 'Gamma prediction', zorder=2)


        ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )


        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)

        if (row_idx == 0) :
            ax.set_ylabel('SNV prevalence', fontsize=16)


        if column_idx == len(nested_species_list)-1:
            ax.set_xlabel('Maximum frequency\nacross hosts, ' + r'$f_{max}$', fontsize=16)


        if column_idx == len(nested_species_list)-2:

            if row_idx > len(nested_species_list[-1])-1:

                ax.set_xlabel('Maximum frequency\nacross hosts, ' + r'$f_{max}$', fontsize=16)

        if (column_idx == 0) and (row_idx == 0):

            ax.legend(loc='upper left', fontsize=14)




fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
fig.savefig("%sf_max_vs_prevalence.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4)
plt.close()
