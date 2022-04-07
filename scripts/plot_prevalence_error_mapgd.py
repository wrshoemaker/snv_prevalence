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

clade_type = 'all'
pi_type = 'pi_include_boundary'


prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()
species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()
n_species_row=5
nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]


def make_plot(variant_type):

    gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
    fig = plt.figure(figsize = (40, 20))


    errors_greter_10_percent_all = []

    for column_idx, column in enumerate(nested_species_list):

        for row_idx, row in enumerate(column):

            ax = fig.add_subplot(gs[column_idx, row_idx])

            predicted_prevalence_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm']
            predicted_prevalence_slm = numpy.asarray(predicted_prevalence_slm)

            observed_prevalence_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
            observed_prevalence_slm = numpy.asarray(observed_prevalence_slm)

            predicted_prevalence_slm_no_zeros = predicted_prevalence_slm[(observed_prevalence_slm>0) & (predicted_prevalence_slm>0) ]
            observed_prevalence_slm_no_zeros = observed_prevalence_slm[(observed_prevalence_slm>0) & (predicted_prevalence_slm>0) ]

            if len(observed_prevalence_slm_no_zeros) == 0:
                continue

            predicted_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd']
            predicted_prevalence = numpy.asarray(predicted_prevalence)

            observed_prevalence = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
            observed_prevalence = numpy.asarray(observed_prevalence)

            predicted_prevalence_no_zeros = predicted_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]
            observed_prevalence_no_zeros = observed_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]


            error_slm = numpy.absolute(predicted_prevalence_slm_no_zeros - observed_prevalence_slm_no_zeros) / observed_prevalence_slm_no_zeros
            error = numpy.absolute(predicted_prevalence_no_zeros - observed_prevalence_no_zeros) / observed_prevalence_no_zeros


            # fraction of SNVs where estimator improves prediction
            error_slm_log10 = numpy.log10(error_slm)
            error_log10 = numpy.log10(error)

            error_slm_log10_no_nan = error_slm_log10[numpy.isfinite(error_slm_log10)]
            error_log10_no_nan = error_log10[numpy.isfinite(error_log10)]

            error_slm_no_nan = 10**error_slm_log10_no_nan
            error_no_nan = 10**error_log10_no_nan

            errors_greter_10_percent_all.append(sum(error_slm_no_nan>=0.1)/len(error_slm_no_nan))

            min_x = min([min(error_slm_no_nan), min(error_no_nan)])
            max_x = max([max(error_slm_no_nan), max(error_no_nan)])

            x_range = numpy.logspace(numpy.log10(min_x), numpy.log10(max_x), num=100, endpoint=True, base=10.0)

            #print(error_no_nan)

            survival_error = [sum(error_no_nan >= i)/len(error_no_nan) for i in x_range]
            survival_slm_error = [sum(error_slm_no_nan >= i)/len(error_slm_no_nan) for i in x_range]

            survival_error = numpy.asarray(survival_error)
            survival_slm_error = numpy.asarray(survival_slm_error)

            ax.plot(x_range, survival_slm_error, ls='-', lw=3, c=prevalence_utils.variant_color_dict[variant_type], zorder=2, label='SLM')
            #ax.plot(x_range, survival_error, ls='--', lw=3, c=prevalence_utils.variant_color_dict[variant_type], zorder=2, label='Single-locus evolution')


            ax.set_xscale('log', basex=10)
            #ax.set_yscale('log', basey=10)

            ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )

            if (row_idx == 0) and (column_idx==0):
                ax.legend(loc="lower left", fontsize=8)

            if (row_idx == 0):
                ax.set_ylabel('Fraction of sites ' + r'$\geq \epsilon$', fontsize=12)

            if column_idx == len(nested_species_list)-1:
                ax.set_xlabel('Relative error of prevalence prediction, ' + r'$\epsilon$', fontsize=12)

            if column_idx == len(nested_species_list)-2:
                if row_idx > len(nested_species_list[-1])-1:
                    ax.set_xlabel('Relative error of of prevalence prediction, ' + r'$\epsilon$', fontsize=12)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%sprevalence_error_mapgd_%s.pdf" % (config.analysis_directory, variant_type), format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
    plt.close()





if __name__=='__main__':

    for variant_type in ['4D', '1D']:

        print(variant_type)

        make_plot(variant_type)
