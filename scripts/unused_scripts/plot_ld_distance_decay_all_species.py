import sys, os
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria

import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#import calculate_snv_distances
import figure_utils
from math import log10,ceil
#import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.random import randint, multinomial
import parse_HMP_data

import figure_utils
import ld_utils

import matplotlib as mpl


bin_width_exponent = 0.3
variant_types = ['4D','1D']



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

ld_dict_all_species = calculate_linkage_disequilibria.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)

ld_species = list(ld_dict_all_species.keys())

n_cols = 6
species_nested_list  = list(ld_utils.chunks(ld_species, n_cols))
n_rows = len(species_nested_list)


fig = plt.figure(figsize = (10, 8))
gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols)

color_dict = {'4D': '#87CEEB', '1D':'#FF6347'}

for row_i_idx, row_i  in enumerate(species_nested_list):

    for column_j_idx, column_j in enumerate(row_i):


        ld_dict = ld_dict_all_species[column_j]

        ax_i_j = fig.add_subplot(gs[ row_i_idx, column_j_idx])
        ax_i_j.xaxis.set_tick_params(labelsize=6)
        ax_i_j.yaxis.set_tick_params(labelsize=6)
        ax_i_j.set_ylim([0.01 , 1.03])

        ax_i_j.set_title(figure_utils.get_pretty_species_name(column_j, include_number=False), fontsize=6)

        if ( row_i_idx == len(species_nested_list)-1) :
            ax_i_j.set_xlabel('Distance between SNVs, $\ell$',  fontsize = 6)


        if row_i_idx == len(species_nested_list)-2:
            if column_j_idx >= len(species_nested_list[row_i_idx + 1]):
                ax_i_j.set_xlabel('Distance between SNVs, $\ell$',  fontsize = 6)


        if column_j_idx == 0:

            ax_i_j.set_ylabel('Linkage disequilibrium, $\sigma^2_d$',  fontsize = 6)


        for variant_type in variant_types:

            distances = ld_dict[variant_type]['distances']
            rsquareds = ld_dict[variant_type]['rsquareds']

            control_rsquared = ld_dict[variant_type]['control_rsquared']
            color = color_dict[variant_type]

            ax_i_j.loglog(distances, rsquareds,'-',color=color, alpha=0.8)

            ax_i_j.loglog([ld_dict[variant_type]['distances'][-1],6e03], [rsquareds[-1], control_rsquared],':',color=color,zorder=21)

            ax_i_j.loglog([6e03], [control_rsquared],'o',color=color,markersize=3,markeredgewidth=0,zorder=21)

            #ax.loglog(ld_dict_variant['all_distances'], ld_dict_variant['all_rsquareds'],'-',color='0.7',label='All samples, %s' % variant_type)


            #ax.set_xlabel('Distance between SNVs, $\ell$')

            #rsquared_ratio_final = ld_dict['1D']['rsquareds'][-1] / ld_dict['4D']['rsquareds'][-1]
            #control_rsquared_ratio = ld_dict['1D']['control_rsquared'] / ld_dict['4D']['control_rsquared']


        #ax_i_j.semilogy([species_idx,species_idx], [rsquareds[idx_900], rsquareds[idx_9]],'-', color=color)
        #ax_i_j.semilogy([species_idx], [rsquareds[idx_9]],'_',markersize=3,color=color) #,markeredgewidth=0)
        #ax_i_j.semilogy([species_idx],[rsquareds[idx_90]],'_',markersize=3, color=color) #,markeredgewidth=0)
        #ax_i_j.semilogy([species_idx], [rsquareds[idx_900]],'_',markersize=3,color=color) #,markeredgewidth=0)
        #line, = species_axis.semilogy([species_idx,species_idx], [control_rsquared, rsquareds[idx_900]],':',color=color)
        #line.set_dashes((0.5,0.75))
        #ax_i_j.semilogy([species_idx], [control_rsquared],'o',markersize=2,markeredgewidth=0,color=color)

# plot of the LD distance decay relationships for all species


fig.subplots_adjust(hspace=0.4, wspace=0.4)
fig.savefig("%s%s.png" % (config.analysis_directory, 'ld_distance_decay_all_species'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
