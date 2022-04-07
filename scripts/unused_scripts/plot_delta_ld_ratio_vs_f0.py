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


n_permute = 100


ld_dict_file = '%sld_all_species_dict.pickle' % (parse_midas_data.data_directory)



with open(ld_dict_file, 'rb') as handle:
    ld_dict_all_species = pickle.load(handle)




# ignore last f0


#Bacteroides_vulgatus_57955

# spcies = {'std': ,  'delta_ell'}
std_delta_dict = {}

all_species_ld = ld_dict_all_species.keys()

#all_species_ld = all_species_ld[:7]

print(all_species_ld)

#all_species_ld = ['Bacteroides_vulgatus_57955', 'Bacteroides_uniformis_57318', 'Odoribacter_splanchnicus_62174']

for species_name in all_species_ld:

    species_ld_dict = ld_dict_all_species[species_name]

    distances = species_ld_dict.keys()

    std_delta_dict[species_name] = {}
    std_delta_dict[species_name]['distances'] = []
    std_delta_dict[species_name]['mean_delta'] = []
    std_delta_dict[species_name]['mean_delta_all'] = []
    std_delta_dict[species_name]['std_mean_delta'] = []
    std_delta_dict[species_name]['std_mean_delta_all'] = []

    for distance in distances:

        if distance <= 20:
            continue

        f0 = species_ld_dict[distance]['f0']
        ratio_1D_4D = species_ld_dict[distance]['1D_4D_ratio']
        ratio_1D_4D_all = species_ld_dict[distance]['1D_4D_ratio_all']

        f0 = numpy.asarray(f0)
        ratio_1D_4D = numpy.asarray(ratio_1D_4D)
        ratio_1D_4D_all = numpy.asarray(ratio_1D_4D_all)

        f0 = f0[1:]
        ratio_1D_4D = ratio_1D_4D[1:]
        ratio_1D_4D_all = ratio_1D_4D_all[1:]

        delta_f0 = f0[1:] - f0[:-1]
        delta_ratio_1D_4D = ratio_1D_4D[1:] - ratio_1D_4D[:-1]
        delta_ratio_1D_4D_all = ratio_1D_4D_all[1:] - ratio_1D_4D_all[:-1]

        delta_ratio_1D_4D_per_f0 = delta_ratio_1D_4D / delta_f0
        delta_ratio_1D_4D_per_f0_all = delta_ratio_1D_4D_all / delta_f0

        mean_delta_ratio_1D_4D_per_f0 = numpy.mean(delta_ratio_1D_4D_per_f0)
        mean_delta_ratio_1D_4D_per_f0_all = numpy.mean(delta_ratio_1D_4D_per_f0_all)

        null_delta_ratio = []
        null_delta_ratio_all = []
        for i in range(n_permute):

            numpy.random.shuffle(ratio_1D_4D)
            numpy.random.shuffle(ratio_1D_4D_all)

            delta_ratio_1D_4D_null = ratio_1D_4D[1:] - ratio_1D_4D[:-1]
            delta_ratio_1D_4D_null_all = ratio_1D_4D_all[1:] - ratio_1D_4D_all[:-1]

            delta_ratio_1D_4D_per_f0_null = delta_ratio_1D_4D_null / delta_f0
            delta_ratio_1D_4D_per_f0_null_all = delta_ratio_1D_4D_null_all / delta_f0

            delta_ratio_1D_4D_per_f0_null_mean = numpy.mean(delta_ratio_1D_4D_per_f0_null)
            delta_ratio_1D_4D_per_f0_null_all_mean = numpy.mean(delta_ratio_1D_4D_per_f0_null_all)

            null_delta_ratio.append(delta_ratio_1D_4D_per_f0_null_mean)
            null_delta_ratio_all.append(delta_ratio_1D_4D_per_f0_null_all_mean)

        null_delta_ratio = numpy.asarray(null_delta_ratio)
        null_delta_ratio_all = numpy.asarray(null_delta_ratio_all)


        std_delta_ratio_1D_4D_per_f0 = (mean_delta_ratio_1D_4D_per_f0 - numpy.mean(null_delta_ratio)) / numpy.std(null_delta_ratio)
        std_delta_ratio_1D_4D_per_f0_all = (mean_delta_ratio_1D_4D_per_f0_all - numpy.mean(null_delta_ratio_all)) / numpy.std(null_delta_ratio_all)


        std_delta_dict[species_name]['distances'].append(distance)


        std_delta_dict[species_name]['mean_delta'].append(mean_delta_ratio_1D_4D_per_f0)
        std_delta_dict[species_name]['mean_delta_all'].append(mean_delta_ratio_1D_4D_per_f0_all)

        std_delta_dict[species_name]['std_mean_delta'].append(std_delta_ratio_1D_4D_per_f0)
        std_delta_dict[species_name]['std_mean_delta_all'].append(std_delta_ratio_1D_4D_per_f0_all)




fig = plt.figure(figsize = (18, 8))

ax_all = plt.subplot2grid((1, 2), (0, 0))
ax_largest_clade = plt.subplot2grid((1, 2), (0, 1))

for species_name in std_delta_dict.keys():

    distances = std_delta_dict[species_name]['distances']
    std_delta = std_delta_dict[species_name]['mean_delta']
    std_delta_all = std_delta_dict[species_name]['mean_delta_all']

    ax_all.scatter(distances, std_delta_all, alpha= 0.3, s=4)
    ax_largest_clade.scatter(distances, std_delta, alpha= 0.3, s=4)




ax_all.set_ylim([-15, 15])
ax_largest_clade.set_ylim([-15, 15])



ax_all.set_xscale('log', basex=10)
ax_largest_clade.set_xscale('log', basex=10)

#ax.set_yscale('log', basey=10)



ax_all.set_xlabel('Distance, ' + r'$\Delta \ell$', fontsize=12)
ax_all.set_ylabel(r'$\left \langle \frac{\Delta \sigma^{2}_{d, N}/\sigma^{2}_{d, S}}{\Delta f_{0}} \right \rangle$', fontsize=11)

ax_largest_clade.set_xlabel('Distance, ' + r'$\Delta \ell$', fontsize=12)
ax_largest_clade.set_ylabel(r'$\left \langle \frac{\Delta \sigma^{2}_{d, N}/\sigma^{2}_{d, S}}{\Delta f_{0}} \right \rangle$', fontsize=11)



ax_all.tick_params(axis = 'both', which = 'major', labelsize = 10)
ax_all.tick_params(axis = 'both', which = 'minor', labelsize = 10)

ax_largest_clade.tick_params(axis = 'both', which = 'major', labelsize = 10)
ax_largest_clade.tick_params(axis = 'both', which = 'minor', labelsize = 10)


ax_all.axhline(y=0, color='k', linestyle=':', lw = 3, zorder=1)
ax_largest_clade.axhline(y=0, color='k', linestyle=':', lw = 3, zorder=1)


ax_all.set_title(figure_utils.clade_type_label_dict['all'] )
ax_largest_clade.set_title(figure_utils.clade_type_label_dict['largest_clade'] )


fig.subplots_adjust(hspace=0.4, wspace=0.2)
fig.savefig("%sdelta_ld_ratio_vs_f0.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
