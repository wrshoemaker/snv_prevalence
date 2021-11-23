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
import matplotlib as mpl

import calculate_coprevalence_f0
import calculate_coprevalence

import calculate_linkage_disequilibria_f0
import calculate_linkage_disequilibria


variant_types = ['4D','1D']
clade_types = ['all','largest_clade']



bin_width_exponent = 1.1

coprevalence_directory = '%scoprevalence_f0/' % (parse_midas_data.data_directory)

ld_directory = '%slinkage_disequilibria_f0/' % (parse_midas_data.data_directory)


list_files_co = [f for f in os.listdir(coprevalence_directory) if os.path.isfile(os.path.join(coprevalence_directory, f))]
list_files_ld = [f for f in os.listdir(ld_directory) if os.path.isfile(os.path.join(ld_directory, f))]



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



all_species_co = set([ "_".join(x.split("_", 3)[:3]) for x in list_files_co])
all_species_ld = set([ "_".join(x.split("_", 3)[:3]) for x in list_files_ld])

all_species = list(all_species_co & all_species_ld)
all_species.remove('Odoribacter_splanchnicus_62174')
all_species.remove('Bacteroidales_bacterium_58650')



all_species = all_species[13:]


#all_species = ['Bacteroidales_bacterium_58650']



for species_name in all_species:

    print(species_name)

    species_co_f0 = [float('0.%s' % re.split('[_ .]+', x)[5]) for x in list_files_co if species_name in x]
    species_co_f0.sort()

    species_ld_f0 = [float('0.%s' % re.split('[_ .]+', x)[5]) for x in list_files_ld if species_name in x]
    species_ld_f0.sort()

    distance_co_dict = {}
    distance_ld_dict = {}

    co_mean_1D_4D_ratio_list = []
    co_mean_1D_4D_ratio_all_list = []

    ld_mean_1D_4D_ratio_list = []
    ld_mean_1D_4D_ratio_all_list = []

    for species_co_f0_i_idx, species_co_f0_i in enumerate(species_co_f0):

        print(species_co_f0_i)

        co_dict_all_species = calculate_coprevalence_f0.calculate_ld_dict([species_name], subject_sample_map, f0=species_co_f0_i, ld_moment=2, bin_width_exponent=bin_width_exponent)
        co_dict_divergence = {}

        co_dict_divergence['1D'] = co_dict_all_species[species_name]['1D']
        co_dict_divergence['4D'] = co_dict_all_species[species_name]['4D']
        co_filter_distances, co_1D_4D_ratio_filter_rsquareds, co_1D_4D_ratio_filter_all_rsquareds, co_ontrol_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(co_dict_divergence)

        if len(co_filter_distances) == 0:
            continue

        co_mean_1D_4D_ratio_list.append(numpy.mean(co_1D_4D_ratio_filter_rsquareds))
        co_mean_1D_4D_ratio_all_list.append(numpy.mean(co_1D_4D_ratio_filter_all_rsquareds))

        for co_filter_distances_j_idx, co_filter_distances_j in enumerate(co_filter_distances):

            if co_filter_distances_j not in distance_co_dict:
                distance_co_dict[co_filter_distances_j] = {}
                distance_co_dict[co_filter_distances_j]['f0'] = []
                distance_co_dict[co_filter_distances_j]['1D_4D_ratio'] = []
                distance_co_dict[co_filter_distances_j]['1D_4D_ratio_all'] = []

            distance_co_dict[co_filter_distances_j]['f0'].append(species_co_f0_i)
            distance_co_dict[co_filter_distances_j]['1D_4D_ratio'].append(co_1D_4D_ratio_filter_rsquareds[co_filter_distances_j_idx])
            distance_co_dict[co_filter_distances_j]['1D_4D_ratio_all'].append(co_1D_4D_ratio_filter_all_rsquareds[co_filter_distances_j_idx])



    #species_ld_f0 = [species_ld_f0[10]]

    #for species_ld_f0_i_idx, species_ld_f0_i in enumerate(species_ld_f0):

    #    ld_dict_all_species = calculate_linkage_disequilibria_f0.calculate_ld_dict([species_name], subject_sample_map, f0=species_ld_f0_i, ld_moment=2, bin_width_exponent=bin_width_exponent)
    #    ld_dict_divergence = {}

        #print(ld_dict_all_species)




    distances_co = list(distance_co_dict.keys())
    distances_co.sort()


    fig = plt.figure(figsize = (8, 8))


    ax_ld_all = plt.subplot2grid((2, 2), (1, 0), colspan=1)
    ax_ld_largest_clade = plt.subplot2grid((2, 2), (1, 1), colspan=1)
    cmap_offset_co = int(0.2*len(distances_co))
    rgb_blue_co = cm.Blues(numpy.linspace(0,1,len(distances_co)+cmap_offset_co ))
    rgb_blue_co = mpl.colors.ListedColormap(rgb_blue_co[cmap_offset_co:,:-1])

    co_ratio_list = ['1D_4D_ratio_all', '1D_4D_ratio']
    for co_ratio_idx, co_ratio in enumerate(co_ratio_list):

        ax_co_ratio = plt.subplot2grid((2, 2), (0, co_ratio_idx), colspan=1)

        #rgb_blue_range
        for distance_j_idx, distance_j in enumerate(distances_co):

            distance_j_f0 = distance_co_dict[distance_j]['f0']
            distance_j_1D_4D_ratio = distance_co_dict[distance_j][co_ratio]

            distance_j_f0 = numpy.asarray(distance_j_f0)
            distance_j_1D_4D_ratio = numpy.asarray(distance_j_1D_4D_ratio)
            _argsort = distance_j_f0.argsort()

            distance_j_f0 = distance_j_f0[_argsort]
            distance_j_1D_4D_ratio = distance_j_1D_4D_ratio[_argsort]


            #ax.plot(distance_j_f0, distance_j_1D_4D_ratio, ls='--', markersize=2, c=rgb_blue[distance_j_idx], marker='o', label='Distance block', alpha=1, zorder=2)
            ax_co_ratio.plot(distance_j_f0, distance_j_1D_4D_ratio, ls='--', markersize=2, c=rgb_blue_co(distance_j_idx), marker='o', alpha=1, zorder=2)



        ax_co_ratio.set_xscale('log', basex=10)
        ax_co_ratio.set_yscale('log', basey=10)

        ax_co_ratio.set_xlabel('Prevalence weight, ' + r'$p_{0}$', fontsize=10)
        ax_co_ratio.set_ylabel('Nonsynonymous/synonymous coprevalence', fontsize=9)

        ax_co_ratio.axhline(y=1, color='k', linestyle=':', lw = 3, zorder=1, label="Neutral model")

        #ax.xaxis.set_tick_params(labelsize=4)
        #ax.yaxis.set_tick_params(labelsize=4)


        #ax.set_xticklabels(fontsize=8)
        #ax.set_yticklabels(fontsize=8)

        ax_co_ratio.set_xlim([min(species_co_f0), max(species_co_f0)])
        #ax.set_ylim([0.1, 1.1])

        #ax.set_xticklabels(ax.get_xticklabels(), fontsize=5)


        ax_co_ratio.tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax_co_ratio.tick_params(axis = 'both', which = 'minor', labelsize = 6)

        ax_co_ratio.set_title(figure_utils.clade_type_label_dict[clade_types[co_ratio_idx]] )



        # colorbar
        norm = mpl.colors.LogNorm(vmin=min(distances_co),vmax=max(distances_co))
        sm = plt.cm.ScalarMappable(cmap=rgb_blue_co, norm=norm)
        sm.set_array([])
        clb = plt.colorbar(sm)

        clb.ax.set_title(r'$\Delta \ell$')





    #ax.legend(loc="upper right", prop={'size': 10})

    #ax.set_title(species_name)
    #fig.text(0, 1.05, figure_utils.get_pretty_species_name(species_name), va='center', fontweight='bold', fontsize=14)
    plt.suptitle(figure_utils.get_pretty_species_name(species_name), fontweight='bold', fontsize=12)

    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    fig.savefig("%sp0_vs_coprevalence/%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
