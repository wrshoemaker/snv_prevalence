import sys, os
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria

#import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

#import calculate_snv_distances
import figure_utils
from math import log10,ceil
from numpy.random import randint, multinomial
import parse_HMP_data

import figure_utils
import ld_utils
import matplotlib as mpl

#bin_width_exponent = 0.3
variant_types = ['4D','1D']

bin_width_exponent = [0.5, 0.7, 0.9, 1.1]

line_styles = ['-', '--', '-.', ':']


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize




sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species_list = parse_midas_data.parse_good_species_list()

ld_dict_all_species_all_bins = {}

for bin_width_exponent in bin_width_exponents:

    ld_dict_all_species = calculate_linkage_disequilibria.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)

    ld_dict_all_species_all_bins[bin_width_exponent] = ld_dict_all_species


ld_species = list(ld_dict_all_species_all_bins[bin_width_exponents[0]].keys())

n_cols = 6
species_nested_list  = list(ld_utils.chunks(ld_species, n_cols))
n_rows = len(species_nested_list)


fig = plt.figure(figsize = (10, 8))
gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols)

plt.rcParams['xtick.labelsize']=4
plt.rcParams['ytick.labelsize']=4


for row_i_idx, row_i  in enumerate(species_nested_list):

    for column_j_idx, column_j in enumerate(row_i):


        ax_i_j = fig.add_subplot(gs[ row_i_idx, column_j_idx])
        #ax_i_j.xaxis.set_tick_params(labelsize=6)
        #ax_i_j.yaxis.set_tick_params(labelsize=6)
        ax_i_j.axhline(y=1, ls=':', color='grey')
        ax_i_j.set_title(figure_utils.get_pretty_species_name(column_j, include_number=False), fontsize=6)

        if ( row_i_idx == len(species_nested_list)-1) :
            ax_i_j.set_xlabel('Distance between SNVs, $\ell$',  fontsize = 6)


        if row_i_idx == len(species_nested_list)-2:
            if column_j_idx >= len(species_nested_list[row_i_idx + 1]):
                ax_i_j.set_xlabel('Distance between SNVs, $\ell$',  fontsize = 6)


        if column_j_idx == 0:

            ax_i_j.set_ylabel('Linkage disequilibrium\nratio, $\sigma^2_{d,N} / \sigma^2_{d,S}$', fontsize = 6)



        for bin_width_exponent_idx, bin_width_exponent in enumerate(bin_width_exponents):

            ld_dict = ld_dict_all_species_all_bins[bin_width_exponent][column_j]

            ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, control_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(ld_dict)

            #ax_i_j.loglog(ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, line_styles[bin_width_exponent_idx], color='k', alpha=0.8)
            #ax_i_j.loglog([  ld_filter_distances[-1],  6e03], [ld_1D_4D_ratio_filter_rsquareds[-1], control_rsquared_ratio],':',color='k',zorder=21)
            #ax_i_j.loglog([6e03], [control_rsquared_ratio],'o',color='k',markersize=3,markeredgewidth=0,zorder=21)

            ax_i_j.plot(ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, line_styles[bin_width_exponent_idx], color='k', alpha=0.8)
            ax_i_j.plot([  ld_filter_distances[-1],  6e03], [ld_1D_4D_ratio_filter_rsquareds[-1], control_rsquared_ratio],':',color='k',zorder=21)
            ax_i_j.plot([6e03], [control_rsquared_ratio],'o',color='k',markersize=3,markeredgewidth=0,zorder=21)

            if bin_width_exponent_idx == 0:
                ax_i_j.set_ylim([min(ld_1D_4D_ratio_filter_rsquareds)*0.9, max(ld_1D_4D_ratio_filter_rsquareds)*1.3])


            # dont know why Line2D isnt working so lets try this
            if (row_i_idx==0) and (column_j_idx==0):
                #ax_i_j.plot([5000], [1000], line_styles[bin_width_exponent_idx], color='k', label= r'$\Delta \ell= 10^{{{}}}$' + str(round(10**(bin_width_exponent),2)) )
                ax_i_j.plot([5000], [1000], line_styles[bin_width_exponent_idx], color='k', label= r'$\Delta \ell= 10^{{{}}}$'.format(str(bin_width_exponent)))

                #ax.text(0.3,0.9, r'$\sigma^{2}_{x} = {{{}}} * x^{{{}}}$'.format(str(round(10**intercept, 3)),  str(round(slope, 3)) ), fontsize=11, color='k', ha='center', va='center', transform=ax.transAxes  )


        if (row_i_idx==0) and (column_j_idx==0):
            ax_i_j.legend(loc='upper left',frameon=False, fontsize=5.5)


        ax_i_j.set_xlim([1, 3006])

        ax_i_j.xaxis.set_tick_params(labelsize=4)
        ax_i_j.yaxis.set_tick_params(labelsize=4)

        ax_i_j.set_xscale('log', basex=10)
        #ax_i_j.set_yscale('log', basey=10)



        #legend_elements_hyper = [Line2D([0], [0], ls='-', color='k', lw=1.5, label=  "444" ),
        #Line2D([0], [0], ls='--', color='k', lw=1.5, label=  r'$\Delta \ell=$' + str(round(10**(bin_width_exponent),3)) ),
        #Line2D([0], [0], ls='-.', color='k', lw=1.5, label=  r'$\Delta \ell=$' + str(round(10**(0.9),3)) ),
        #Line2D([0], [0], ls=':', color='k', lw=1.5, label=  r'$\Delta \ell=$' + str(round(10**(1.1),3)) )]

        #ax_i_j.legend(legend_elements_hyper, loc='upper right', fontsize=7.5)




plt.rcParams['xtick.labelsize']=4
plt.rcParams['ytick.labelsize']=4





fig.subplots_adjust(hspace=0.4, wspace=0.4)
fig.savefig("%s%s.png" % (config.analysis_directory, 'ld_ratio_distance_decay_all_species'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
