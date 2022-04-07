import sys, os
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria_divergence

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
from matplotlib import cm

#bin_width_exponent = 0.3
variant_types = ['4D','1D']

bin_width_exponent = 1.1



#color_range =  numpy.linspace(0.0, 1.0, 18)
#rgb_blue = cm.get_cmap('Blues')( color_range )


line_styles = ['-', '--', '-.', ':', ':']


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
#good_species_list = ['Bacteroides_vulgatus_57955']
ld_dict_all_species = calculate_linkage_disequilibria_divergence.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)
#print(ld_dict_all_species.keys())

ld_species = list(ld_dict_all_species.keys())

n_cols = 6
species_nested_list  = list(ld_utils.chunks(ld_species, n_cols))
n_rows = len(species_nested_list)



fig = plt.figure(figsize = (10, 8))
gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols)

plt.rcParams['xtick.labelsize']=4
plt.rcParams['ytick.labelsize']=4



#divergences = []

for row_i_idx, row_i  in enumerate(species_nested_list):

    for column_j_idx, column_j in enumerate(row_i):

        ld_dict = ld_dict_all_species[column_j]

        divergences = list(set([x[0] for x in list(ld_dict.keys()) if 'divergence' in x[0]]))
        divergences_float = [float(x.split('_')[1]) for x in divergences]

        # sort the strings by float value
        zip_divergences = zip(divergences, divergences_float)
        zip_divergences = sorted(zip_divergences, key=lambda x: x[1])

        ax_i_j = fig.add_subplot(gs[ row_i_idx, column_j_idx])

        ax_i_j.axhline(y=1, ls=':', color='grey')
        ax_i_j.set_title(figure_utils.get_pretty_species_name(column_j, include_number=False), fontsize=6)

        if ( row_i_idx == len(species_nested_list)-1) :
            ax_i_j.set_xlabel('Divergence, $d$',  fontsize = 6)


        if row_i_idx == len(species_nested_list)-2:
            if column_j_idx >= len(species_nested_list[row_i_idx + 1]):
                ax_i_j.set_xlabel('Divergence, $d$',  fontsize = 6)


        if column_j_idx == 0:
            ax_i_j.set_ylabel('Linkage disequilibrium\nratio, $\sigma^2_{d,N} / \sigma^2_{d,S}$', fontsize = 6)


        divergence_vs_ratio_dict = {}
        divergence_float_all = []
        for zip_divergence_idx, zip_divergence in enumerate(zip_divergences):

            divergence = zip_divergence[0]
            divergence_float = zip_divergence[1]
            divergence_float_all.append(divergence_float)

            #divergence_float_list.append(divergence_float)

            ld_dict_divergence_1D = ld_dict[(divergence, '1D')]
            ld_dict_divergence_4D = ld_dict[(divergence, '4D')]

            ld_dict_divergence = {}
            ld_dict_divergence['1D'] = ld_dict_divergence_1D
            ld_dict_divergence['4D'] = ld_dict_divergence_4D

            ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, control_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(ld_dict_divergence)

            for ld_filter_distance_idx, ld_filter_distance in enumerate(ld_filter_distances):

                if ld_filter_distance not in divergence_vs_ratio_dict:
                    divergence_vs_ratio_dict[ld_filter_distance] = {}
                    divergence_vs_ratio_dict[ld_filter_distance]['ld_1D_4D_ratio'] = []
                    divergence_vs_ratio_dict[ld_filter_distance]['divergence_float'] = []

                divergence_vs_ratio_dict[ld_filter_distance]['ld_1D_4D_ratio'].append(ld_1D_4D_ratio_filter_rsquareds[ld_filter_distance_idx])
                divergence_vs_ratio_dict[ld_filter_distance]['divergence_float'].append(divergence_float)


        for distance_chunk in divergence_vs_ratio_dict.keys():
            if len(divergence_vs_ratio_dict[distance_chunk]['ld_1D_4D_ratio']) != len(divergences):
                del divergence_vs_ratio_dict[distance_chunk]


        all_ratios = [ divergence_vs_ratio_dict[distance]['ld_1D_4D_ratio'] for distance in divergence_vs_ratio_dict.keys()]
        all_ratios = numpy.array(all_ratios)

        all_ratios_mean = numpy.mean(all_ratios, axis=0)



        for distance_chunk in divergence_vs_ratio_dict.keys():

            distance_chunk_x = divergence_vs_ratio_dict[distance_chunk]['divergence_float']
            distance_chunk_y = divergence_vs_ratio_dict[distance_chunk]['ld_1D_4D_ratio']

            ax_i_j.plot(distance_chunk_x, distance_chunk_y, '-', color='grey', alpha=0.05, zorder=1)


        ax_i_j.plot(divergence_float_all, all_ratios_mean, ls='--', marker='o', color='k', alpha=1, zorder=2)

        #ax_i_j.set_xlim([1, 3006])

        ax_i_j.xaxis.set_tick_params(labelsize=4)
        ax_i_j.yaxis.set_tick_params(labelsize=4)

        ax_i_j.set_xscale('log', basex=10)



fig.subplots_adjust(hspace=0.4, wspace=0.4)
fig.savefig("%s%s.png" % (config.analysis_directory, 'ld_ratio_chunk_vs_divergence'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
