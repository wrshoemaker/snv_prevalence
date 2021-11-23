from __future__ import division
import sample_utils
import config
import parse_midas_data
import os.path
import os
import re
import pylab
import sys
import numpy
import gzip
import pickle
import bz2

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
import ld_utils
from math import log10,ceil,fabs
from numpy.random import randint, choice, multinomial

import midas_db_utils

import scipy.stats as stats

import parse_HMP_data
import figure_utils

#good_species_list = [good_species_list[3]]

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import calculate_coprevalence_f0
import calculate_coprevalence


variant_types = ['4D','1D']

bin_width_exponent = 1.1

ld_directory = '%scoprevalence_f0/' % (parse_midas_data.data_directory)

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
#ld_dict_all_species = calculate_coprevalence.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)
#print(ld_dict_all_species.keys())

#print(list_files)



#all_species = ['Bacteroides_vulgatus_57955']

all_species = list(set([ "_".join(x.split("_", 3)[:3]) for x in list_files]))

fig, ax = plt.subplots(figsize=(4,4))

for species_name in all_species:

    species_f0 = [float('0.%s' % re.split('[_ .]+', x)[5]) for x in list_files if species_name in x]
    species_f0.sort()

    distance_ld_dict = {}

    mean_1D_4D_ratio_list = []

    for species_f0_i_idx, species_f0_i in enumerate(species_f0):


        ld_dict_all_species = calculate_coprevalence_f0.calculate_ld_dict([species_name], subject_sample_map, f0=species_f0_i, ld_moment=2, bin_width_exponent=bin_width_exponent)

        ld_dict_divergence = {}
        ld_dict_divergence['1D'] = ld_dict_all_species[species_name]['1D']
        ld_dict_divergence['4D'] = ld_dict_all_species[species_name]['4D']
        ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, control_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(ld_dict_divergence)

        mean_1D_4D_ratio_list.append(numpy.mean(ld_1D_4D_ratio_filter_rsquareds))

        for ld_filter_distances_j_idx, ld_filter_distances_j in enumerate(ld_filter_distances):

            if ld_filter_distances_j not in distance_ld_dict:
                distance_ld_dict[ld_filter_distances_j] = {}
                distance_ld_dict[ld_filter_distances_j]['f0'] = []
                distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio'] = []

            distance_ld_dict[ld_filter_distances_j]['f0'].append(species_f0_i)
            distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio'].append(ld_1D_4D_ratio_filter_rsquareds[ld_filter_distances_j_idx])




    dict_1D_4D_ratio = {}
    for species_f0_i_idx, species_f0_i in enumerate(species_f0):
        dict_1D_4D_ratio[species_f0_i] = []

    #color_range =  numpy.linspace(0.0, 1.0, len(distances))
    #rgb_blue = cm.get_cmap('Blues')( color_range )
    for distance_j in distance_ld_dict.keys():

        distance_j_f0 = distance_ld_dict[distance_j]['f0']
        distance_j_1D_4D_ratio = distance_ld_dict[distance_j]['1D_4D_ratio']

        distance_j_f0 = numpy.asarray(distance_j_f0)
        distance_j_1D_4D_ratio = numpy.asarray(distance_j_1D_4D_ratio)
        _argsort = distance_j_f0.argsort()

        distance_j_f0 = distance_j_f0[_argsort]
        distance_j_1D_4D_ratio = distance_j_1D_4D_ratio[_argsort]

        for f0_1D_4D_ratio in zip(distance_j_f0, distance_j_1D_4D_ratio):

            dict_1D_4D_ratio[f0_1D_4D_ratio[0]].append(f0_1D_4D_ratio[1])


    print(dict_1D_4D_ratio.keys())

    

    ratio_to_plot = [10**numpy.mean(numpy.log10(dict_1D_4D_ratio[f0])) for f0 in species_f0]


    ax.scatter(distance_j_f0, ratio_to_plot, alpha = 0.8)

    #ax.plot(distance_j_f0, distance_j_1D_4D_ratio, ls='--', markersize=2, c=rgb_blue[distance_j_idx], marker='o', label='Distance block', alpha=0.8, zorder=2)


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel('Prevalence weight', fontsize=12)
ax.set_ylabel('Nonsynonymous/synonymous coprevalence', fontsize=12)

ax.axhline(y=1, color='k', linestyle=':', lw = 3, zorder=1, label="Neutral model")

ax.legend(loc="upper right")

#ax.set_title(species_name)

fig.subplots_adjust(hspace=0.4, wspace=0.4)
fig.savefig("%s%s.png" % (config.analysis_directory, 'p0_vs_coprevalence_ratio'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()






    # old code for all species.
    # get mean ratio across distance bins

    #dict_1D_4D_ratio = {}
    #for species_f0_i_idx, species_f0_i in enumerate(species_f0):
    #    dict_1D_4D_ratio[species_f0_i] = []

    #for distance_j in distance_ld_dict.keys():

        #distance_j_f0 = distance_ld_dict[distance_j]['f0']
        #distance_j_1D_4D_ratio = distance_ld_dict[distance_j]['1D_4D_ratio']

        #distance_j_f0 = numpy.asarray(distance_j_f0)
        #distance_j_1D_4D_ratio = numpy.asarray(distance_j_1D_4D_ratio)
        #_argsort = distance_j_f0.argsort()

        #for f0_1D_4D_ratio in zip(distance_j_f0, distance_j_1D_4D_ratio):

        #    dict_1D_4D_ratio[f0_1D_4D_ratio[0]].append(f0_1D_4D_ratio[1])


        #distance_j_f0 = distance_j_f0[_argsort[::-1]]
        #distance_j_1D_4D_ratio = distance_j_1D_4D_ratio[_argsort[::-1]]

        #distance_j_f0 = distance_j_f0[_argsort]
        #distance_j_1D_4D_ratio = distance_j_1D_4D_ratio[_argsort]


        #zipped_lists = zip(list1, list2)

    #ratio_to_plot = [numpy.mean(dict_1D_4D_ratio[f0]) for f0 in species_f0]


    #ax.scatter(distance_j_f0, distance_j_1D_4D_ratio, alpha = 0.8)

    #ax.plot(species_f0, ratio_to_plot, ls='--', markersize=2, marker='o', alpha=0.3, zorder=2)

    #ax_i_j.plot(distance_j_f0, distance_j_1D_4D_ratio, '-', color='grey', alpha=0.05, zorder=1)

    #ax_i_j.xaxis.set_tick_params(labelsize=4)
    #ax_i_j.yaxis.set_tick_params(labelsize=4)
