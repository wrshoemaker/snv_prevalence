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

#import calculate_coprevalence_f0
#import calculate_coprevalence

import calculate_linkage_disequilibria_f0
import calculate_linkage_disequilibria


variant_types = ['4D','1D']
clade_types = ['all','largest_clade']



bin_width_exponent = 1.1

ld_directory = '%slinkage_disequilibria_f0/' % (parse_midas_data.data_directory)


sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()




list_files_ld = [f for f in os.listdir(ld_directory) if os.path.isfile(os.path.join(ld_directory, f))]

all_species_ld = list(set([ "_".join(x.split("_", 3)[:3]) for x in list_files_ld]))

all_species_ld.remove('Phascolarctobacterium_sp_59817')
all_species_ld.remove('Akkermansia_muciniphila_55290')


#all_species_ld = ['Bacteroides_vulgatus_57955']

#species_name = 'Bacteroides_vulgatus_57955'


ld_dict_all_species_to_pickle = {}



for species_name in all_species_ld:

    print(species_name)


    distance_ld_dict = {}

    ld_mean_1D_4D_ratio_list = []
    ld_mean_1D_4D_ratio_all_list = []

    species_ld_f0 = [float('0.%s' % re.split('[_ .]+', x)[5]) for x in list_files_ld if species_name in x]
    species_ld_f0.sort()


    #species_ld_f0 = [species_ld_f0[-1]]

    for species_ld_f0_i_idx, species_ld_f0_i in enumerate(species_ld_f0):


        ld_dict_all_species = calculate_linkage_disequilibria_f0.calculate_ld_dict([species_name], subject_sample_map, f0=species_ld_f0_i, ld_moment=2, bin_width_exponent=bin_width_exponent)

        ld_dict_divergence = {}

        if len(ld_dict_all_species[species_name]['1D']) != len(ld_dict_all_species[species_name]['4D']):
            continue

        print(species_ld_f0_i)

        ld_dict_divergence['1D'] = ld_dict_all_species[species_name]['1D']
        ld_dict_divergence['4D'] = ld_dict_all_species[species_name]['4D']
        ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, ld_1D_4D_ratio_filter_all_rsquareds, ld_control_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(ld_dict_divergence)

        if len(ld_filter_distances) == 0:
            continue

        ld_mean_1D_4D_ratio_list.append(numpy.mean(ld_1D_4D_ratio_filter_rsquareds))
        ld_mean_1D_4D_ratio_all_list.append(numpy.mean(ld_1D_4D_ratio_filter_all_rsquareds))

        for ld_filter_distances_j_idx, ld_filter_distances_j in enumerate(ld_filter_distances):

            if ld_filter_distances_j not in distance_ld_dict:
                distance_ld_dict[ld_filter_distances_j] = {}
                distance_ld_dict[ld_filter_distances_j]['f0'] = []
                distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio'] = []
                distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio_all'] = []

            distance_ld_dict[ld_filter_distances_j]['f0'].append(species_ld_f0_i)
            distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio'].append(ld_1D_4D_ratio_filter_rsquareds[ld_filter_distances_j_idx])
            distance_ld_dict[ld_filter_distances_j]['1D_4D_ratio_all'].append(ld_1D_4D_ratio_filter_all_rsquareds[ld_filter_distances_j_idx])



    distances_ld = list(distance_ld_dict.keys())
    distances_ld.sort()


    if len(distance_ld_dict) == 0:
        continue


    ld_dict_all_species_to_pickle[species_name] = distance_ld_dict



    fig = plt.figure(figsize = (18, 8))

    #ax_ld_all = plt.subplot2grid((2, 2), (1, 0), colspan=1)
    #ax_ld_largest_clade = plt.subplot2grid((2, 2), (1, 1), colspan=1)
    cmap_offset_ld = int(0.2*len(distances_ld))
    rgb_blue_ld = cm.Blues(numpy.linspace(0,1,len(distances_ld)+cmap_offset_ld ))
    rgb_blue_ld = mpl.colors.ListedColormap(rgb_blue_ld[cmap_offset_ld:,:-1])

    ld_ratio_list = ['1D_4D_ratio_all', '1D_4D_ratio']
    for ld_ratio_idx, ld_ratio in enumerate(ld_ratio_list):

        ax_ld_ratio = plt.subplot2grid((1, 2), (0, ld_ratio_idx), colspan=1)

        #rgb_blue_range
        f0_all = []
        for distance_j_idx, distance_j in enumerate(distances_ld):

            distance_j_f0 = distance_ld_dict[distance_j]['f0']
            distance_j_1D_4D_ratio = distance_ld_dict[distance_j][ld_ratio]

            distance_j_f0 = numpy.asarray(distance_j_f0)
            distance_j_1D_4D_ratio = numpy.asarray(distance_j_1D_4D_ratio)
            _argsort = distance_j_f0.argsort()

            distance_j_f0 = distance_j_f0[_argsort]
            distance_j_1D_4D_ratio = distance_j_1D_4D_ratio[_argsort]

            f0_all.extend(distance_j_f0.tolist())


            ax_ld_ratio.plot(distance_j_f0, distance_j_1D_4D_ratio, ls='--', markersize=2, c=rgb_blue_ld(distance_j_idx), marker='o', alpha=1, zorder=2)



        ax_ld_ratio.set_xscale('log', basex=10)
        ax_ld_ratio.set_yscale('log', basey=10)

        ax_ld_ratio.set_xlabel('Frequency weight, ' + r'$f_{0}$', fontsize=10)
        ax_ld_ratio.set_ylabel('Nonsynonymous/synonymous LD', fontsize=9)

        ax_ld_ratio.axhline(y=1, color='k', linestyle=':', lw = 3, zorder=1, label="Neutral model")



        ax_ld_ratio.set_xlim([min(f0_all), max(f0_all)])
        #ax_ld_ratio.set_xlim([0.03, max(species_ld_f0)])

        if ld_ratio == '1D_4D_ratio_all':

            ax_ld_ratio.set_ylim([0.5, 1.3])

        #else:

        #    ax_ld_ratio.set_ylim([0.5, 1.3])


        ax_ld_ratio.tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax_ld_ratio.tick_params(axis = 'both', which = 'minor', labelsize = 6)

        ax_ld_ratio.set_title(figure_utils.clade_type_label_dict[clade_types[ld_ratio_idx]] )



        # colorbar
        norm = mpl.colors.LogNorm(vmin=min(distances_ld),vmax=max(distances_ld))
        sm = plt.cm.ScalarMappable(cmap=rgb_blue_ld, norm=norm)
        sm.set_array([])
        clb = plt.colorbar(sm)

        clb.ax.set_title(r'$\Delta \ell$')





    #ax.legend(loc="upper right", prop={'size': 10})

    #ax.set_title(species_name)
    #fig.text(0, 1.05, figure_utils.get_pretty_species_name(species_name), va='center', fontweight='bold', fontsize=14)
    plt.suptitle(figure_utils.get_pretty_species_name(species_name), fontweight='bold', fontsize=12)

    fig.subplots_adjust(hspace=0.4, wspace=0.2)
    fig.savefig("%sf0_vs_ld/%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()






ld_dict_file = '%sld_all_species_dict.pickle' % (parse_midas_data.data_directory)


with open(ld_dict_file, 'wb') as outfile:
    pickle.dump(ld_dict_all_species_to_pickle, outfile, protocol=pickle.HIGHEST_PROTOCOL)
