import sys, os, re
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria_divergence
import calculate_linkage_disequilibria_f0

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

ld_directory = '%slinkage_disequilibria_f0/' % (parse_midas_data.data_directory)

list_files = [f for f in os.listdir(ld_directory) if os.path.isfile(os.path.join(ld_directory, f))]


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



for row_i_idx, row_i  in enumerate(species_nested_list):

    for column_j_idx, column_j in enumerate(row_i):

        #if (row_i_idx!=0) or (column_j_idx != 0):
        #    continue

        ax_i_j = fig.add_subplot(gs[ row_i_idx, column_j_idx])

        ax_i_j.axhline(y=1, ls=':', color='k')
        ax_i_j.set_title(figure_utils.get_pretty_species_name(column_j, include_number=False), fontsize=6)

        if ( row_i_idx == len(species_nested_list)-1) :
            ax_i_j.set_xlabel('Minimum frequency, $f_{0}$',  fontsize = 6)

        if row_i_idx == len(species_nested_list)-2:
            if column_j_idx >= len(species_nested_list[row_i_idx + 1]):
                ax_i_j.set_xlabel('Minimum frequency, $f_{0}$',  fontsize = 6)

        if column_j_idx == 0:
            ax_i_j.set_ylabel('Linkage disequilibrium\nratio, $\sigma^2_{d,N} / \sigma^2_{d,S}$', fontsize = 6)


        species_f0 = [float('0.%s' % re.split('[_ .]+', x)[5]) for x in list_files if column_j in x]
        species_f0.sort()

        distance_ld_dict = {}

        mean_1D_4D_ratio_list = []

        for species_f0_i_idx, species_f0_i in enumerate(species_f0):

            ld_dict_all_species = calculate_linkage_disequilibria_f0.calculate_ld_dict([column_j], subject_sample_map, f0=species_f0_i, bin_width_exponent=bin_width_exponent)

            ld_dict_divergence = {}
            ld_dict_divergence['1D'] = ld_dict_all_species[column_j]['1D']
            ld_dict_divergence['4D'] = ld_dict_all_species[column_j]['4D']
            ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, control_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(ld_dict_divergence)

            mean_1D_4D_ratio_list.append(numpy.mean(ld_1D_4D_ratio_filter_rsquareds))

            for ld_filter_distances_j_idx, ld_filter_distances_j in enumerate(ld_filter_distances):

                if ld_filter_distances_j not in distance_ld_dict:
                    distance_ld_dict[ld_filter_distances_j] = {}
                    distance_ld_dict[ld_filter_distances_j]['f0'] = []
                    distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio'] = []

                distance_ld_dict[ld_filter_distances_j]['f0'].append(species_f0_i)
                distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio'].append(ld_1D_4D_ratio_filter_rsquareds[ld_filter_distances_j_idx])

        for distance_j in distance_ld_dict.keys():

            distance_j_f0 = distance_ld_dict[distance_j]['f0']
            distance_j_1D_4D_ratio = distance_ld_dict[distance_j]['1D_4D_ratio']

            distance_j_f0 = numpy.asarray(distance_j_f0)
            distance_j_1D_4D_ratio = numpy.asarray(distance_j_1D_4D_ratio)

            _argsort = distance_j_f0.argsort()

            distance_j_f0 = distance_j_f0[_argsort[::-1]]
            distance_j_1D_4D_ratio = distance_j_1D_4D_ratio[_argsort[::-1]]


            ax_i_j.plot(distance_j_f0, distance_j_1D_4D_ratio, '-', color='grey', alpha=0.05, zorder=1)


        ax_i_j.plot(species_f0, mean_1D_4D_ratio_list, ls='--', markersize=2, marker='o', color='k', alpha=1, zorder=2)

        ax_i_j.xaxis.set_tick_params(labelsize=4)
        ax_i_j.yaxis.set_tick_params(labelsize=4)

        ax_i_j.set_xscale('log', basex=10)



fig.subplots_adjust(hspace=0.4, wspace=0.4)
fig.savefig("%s%s.png" % (config.analysis_directory, 'f0_vs_ld_ratio'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
