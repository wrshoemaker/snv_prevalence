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

#clade_types = ['all','largest_clade', 'no_strains']


#prevalence_dict = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()
prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()
species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()
n_species_row=5
nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]

label_dict = {'4D': 'Synonymous', '1D': 'Nonsynonymous'}


clade_type = 'all'
pi_type = 'pi_include_boundary'
#variant_type = '4D'
best = False

if best == True:
    best_status = '_best'
else:
    best_status = ''

variant_types = ['4D', '1D']

gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
fig = plt.figure(figsize = (40, 20))

for column_idx, column in enumerate(nested_species_list):

    for row_idx, row in enumerate(column):

        ax = fig.add_subplot(gs[column_idx, row_idx])

        for variant_type in variant_types:

            n_non_zero_frequencies = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
            n_set = list(set(n_non_zero_frequencies))
            n_non_zero_frequencies = numpy.asarray(n_non_zero_frequencies)

            predicted_prevalence_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm%s' % best_status]
            predicted_prevalence_slm = numpy.asarray(predicted_prevalence_slm)

            observed_prevalence_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm%s' % best_status]
            observed_prevalence_slm = numpy.asarray(observed_prevalence_slm)

            idx_ = (observed_prevalence_slm>0) & (predicted_prevalence_slm>0)

            predicted_prevalence_slm = predicted_prevalence_slm[idx_]
            observed_prevalence_slm = observed_prevalence_slm[idx_]
            n_non_zero_frequencies = n_non_zero_frequencies[idx_]

            error_slm = numpy.absolute(observed_prevalence_slm - predicted_prevalence_slm) / observed_prevalence_slm

            mean_error_slm_n = [numpy.mean(error_slm[n_non_zero_frequencies == n]) for n in n_set]


            ax.plot(n_set, mean_error_slm_n, color=prevalence_utils.variant_color_dict[variant_type], label=label_dict[variant_type])



        ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )

        #ax.set_yscale('log', basey=10)

        if (row_idx == 0) and (column_idx == 0):
            ax.legend(loc="upper right", fontsize=8)


        if (row_idx == 0):
            ax.set_ylabel('Relative error of SLM', fontsize=12)


        if column_idx == len(nested_species_list)-1:
            ax.set_xlabel('Number of hosts where SNV is present', fontsize=12)


        if column_idx == len(nested_species_list)-2:

            if row_idx > len(nested_species_list[-1])-1:

                ax.set_xlabel('Number of hosts where SNV is present', fontsize=12)




fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# dpi = 600
fig.savefig("%sprevalence_vs_error_slm_%s%s.png" % (config.analysis_directory, clade_type, best_status), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
plt.close()
