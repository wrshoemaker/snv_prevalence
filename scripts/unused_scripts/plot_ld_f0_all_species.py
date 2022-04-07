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

ld_directory = '%slinkage_disequilibria_f0_same_range/' % (parse_midas_data.data_directory)

ld_dict_file = '%sld_all_species_dict_same_range.pickle' % (parse_midas_data.data_directory)


sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()




list_files_ld = [f for f in os.listdir(ld_directory) if os.path.isfile(os.path.join(ld_directory, f))]

all_species_ld = list(set([ "_".join(x.split("_", 3)[:3]) for x in list_files_ld]))

all_species_ld.remove('Phascolarctobacterium_sp_59817')
all_species_ld.remove('Akkermansia_muciniphila_55290')

#all_species_ld = ['Bacteroides_vulgatus_57955', 'Bacteroides_uniformis_57318']

ld_all_species_dict = {}


species_count = 0

for species_name in all_species_ld:


    distance_ld_dict = {}

    ld_mean_1D_4D_ratio_list = []
    ld_mean_1D_4D_ratio_all_list = []

    species_ld_f0 = [float('0.%s' % re.split('[_ .]+', x)[5]) for x in list_files_ld if species_name in x]
    species_ld_f0.sort()

    #species_ld_f0 = [species_ld_f0[-1]]

    #species_ld_f0 = [species_ld_f0[-20], species_ld_f0[-10], species_ld_f0[-1] ]



    for species_ld_f0_i_idx, species_ld_f0_i in enumerate(species_ld_f0):


        ld_dict_all_species = calculate_linkage_disequilibria_f0.calculate_ld_dict([species_name], subject_sample_map, f0=species_ld_f0_i, ld_moment=2, same_f0_range=True, bin_width_exponent=bin_width_exponent)

        ld_dict_divergence = {}

        if len(ld_dict_all_species[species_name]['1D']) != len(ld_dict_all_species[species_name]['4D']):
            continue

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



    #ld_ratio_list = ['1D_4D_ratio_all', '1D_4D_ratio']

    #rgb_blue_range
    f0_all = []
    for distance_j_idx, distance_j in enumerate(distances_ld):

        distance_j_f0 = distance_ld_dict[distance_j]['f0']
        distance_j_1D_4D_ratio = distance_ld_dict[distance_j]['1D_4D_ratio']

        #distance_j_f0 = numpy.asarray(distance_j_f0)
        #distance_j_1D_4D_ratio = numpy.asarray(distance_j_1D_4D_ratio)
        #_argsort = distance_j_f0.argsort()

        #distance_j_f0 = distance_j_f0[_argsort]
        #distance_j_1D_4D_ratio = distance_j_1D_4D_ratio[_argsort]



        if distance_j not in ld_all_species_dict:
            ld_all_species_dict[distance_j] = {}

        #if distance_j_1D_4D_ratio[0] > 1:

        #    ld_all_species_dict[distance_j]['N_greater'] += 1

        #ld_all_species_dict[distance_j]['N'] += 1


        for distance_j_f0_i, distance_j_1D_4D_ratio_i in zip(distance_j_f0, distance_j_1D_4D_ratio):

            if distance_j_f0_i not in ld_all_species_dict[distance_j]:
                ld_all_species_dict[distance_j][distance_j_f0_i] = {}
                ld_all_species_dict[distance_j][distance_j_f0_i]['N'] = 0
                ld_all_species_dict[distance_j][distance_j_f0_i]['N_greater'] = 0

                ld_all_species_dict[distance_j][distance_j_f0_i]['ld_ratio'] = []

                ld_all_species_dict[distance_j][distance_j_f0_i]['species'] = []


            if distance_j_1D_4D_ratio_i > 1:

                ld_all_species_dict[distance_j][distance_j_f0_i]['N_greater'] += 1

            ld_all_species_dict[distance_j][distance_j_f0_i]['ld_ratio'].append(distance_j_1D_4D_ratio_i)

            ld_all_species_dict[distance_j][distance_j_f0_i]['species'].append(species_name)

            ld_all_species_dict[distance_j][distance_j_f0_i]['N'] += 1





with open(ld_dict_file, 'wb') as outfile:
    pickle.dump(ld_all_species_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)



distance_ld_all_species_dict = ld_all_species_dict.keys()

cmap_offset_ld = int(0.2*len(distance_ld_all_species_dict))
rgb_blue_ld = cm.Blues(numpy.linspace(0,1,len(distance_ld_all_species_dict)+cmap_offset_ld ))
rgb_blue_ld = mpl.colors.ListedColormap(rgb_blue_ld[cmap_offset_ld:,:-1])


fig, ax = plt.subplots(figsize=(4,4))

for distance_j_idx, distance_j in enumerate(distance_ld_all_species_dict):

    f0_range = ld_all_species_dict[distance_j].keys()
    f0_range.sort()

    fraction_greater = []

    mean_ld_ratio = []

    for f0_range_i in f0_range:

        fraction_greater_i = ld_all_species_dict[distance_j][f0_range_i]['N_greater'] / ld_all_species_dict[distance_j][f0_range_i]['N']

        fraction_greater.append(fraction_greater_i)

        ld_ratio = ld_all_species_dict[distance_j][f0_range_i]['ld_ratio']


        mean_ld_ratio_i = 10**numpy.mean(numpy.log10(ld_ratio))

        mean_ld_ratio.append(mean_ld_ratio_i)


    ax.plot(f0_range, mean_ld_ratio, ls='--', markersize=2, c=rgb_blue_ld(distance_j_idx), marker='o', alpha=0.8, zorder=2)



ax.set_ylim([10**-0.05, 10**0.05])

ax.axhline(y=1, color='k', linestyle=':', lw = 3, zorder=3)

ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)



ax.set_xlabel('Frequency weight, ' + r'$f_{0}$', fontsize=10)
ax.set_ylabel('Mean ' + r'$\sigma^{2}_{d, N}/\sigma^{2}_{d, S}$' + ' across species' , fontsize=9)



ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
ax.tick_params(axis = 'both', which = 'minor', labelsize = 6)

# colorbar
norm = mpl.colors.LogNorm(vmin=min(distance_ld_all_species_dict),vmax=max(distance_ld_all_species_dict))
sm = plt.cm.ScalarMappable(cmap=rgb_blue_ld, norm=norm)
sm.set_array([])
clb = plt.colorbar(sm)

clb.ax.set_title(r'$\Delta \ell$')

#    ld_all_species_dict[distance_j]['fraction_greater'] = ld_all_species_dict[distance_j]['N_greater'] / ld_all_species_dict[distance_j]['N']


fig.subplots_adjust(hspace=0.4, wspace=0.2)
fig.savefig("%sf0_ld_all_species.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
