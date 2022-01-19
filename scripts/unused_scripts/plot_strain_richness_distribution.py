from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import collections
import os.path
import scipy.stats as stats

import diversity_utils
import core_gene_utils
import parse_midas_data
import parse_HMP_data
import sample_utils
import calculate_substitution_rates
import clade_utils
import plot_utils
import figure_utils

import matplotlib.pyplot as plt

from scipy.stats import gamma
import scipy.special
from scipy.integrate import quad

import calculate_predicted_prevalence_mapgd

good_species_list = calculate_predicted_prevalence_mapgd.good_species_list

species_color_map, ordered_species_list = plot_utils.get_species_color_map()


pi_type = 'pi_include_boundary'
variant_type = '4D'

max_richness = 4
richness_range = list(range(1, max_richness+1))



#fig, ax = plt.subplots(figsize=(4,4))

fig = plt.figure(figsize = (4, 8)) #
fig.subplots_adjust(bottom= 0.15)

ax_samples = plt.subplot2grid((2, 1), (0, 0), colspan=1)
ax_strains = plt.subplot2grid((2, 1), (1, 0), colspan=1)



relative_richness_dict, number_samples_dict = prevalence_utils.get_relative_richness_dict()

# plot number samples
number_samples = [number_samples_dict[s] for s in good_species_list]
good_species_list_sorted_samples = [s[0] for s in sorted(zip(good_species_list, number_samples), key = lambda t: t[1])]

good_species_list_sorted_samples_pretty = [figure_utils.get_pretty_species_name(s) for s in good_species_list_sorted_samples][::-1]
sorted_samples = [s[1] for s in sorted(zip(good_species_list, number_samples), key = lambda t: t[1])][::-1]

ax_samples.barh(good_species_list_sorted_samples_pretty, sorted_samples, height=0.8, align='center', color='dodgerblue')
ax_samples.set_xlabel('Number of hosts', fontsize=11)

ax_samples.xaxis.set_tick_params(labelsize=8)
ax_samples.yaxis.set_tick_params(labelsize=7)
ax_samples.set_ylim([-0.6, len(good_species_list_sorted_samples)-0.3])





good_species_list_sorted_pretty = [figure_utils.get_pretty_species_name(s) for s in good_species_list_sorted]


proprtion_richness_1, proprtion_richness_2, proprtion_richness_3, proprtion_richness_4  = proprtion_richness

proprtion_richness_1 = numpy.asarray(proprtion_richness_1)
proprtion_richness_2 = numpy.asarray(proprtion_richness_2)
proprtion_richness_3 = numpy.asarray(proprtion_richness_3)
proprtion_richness_4 = numpy.asarray(proprtion_richness_4)

ax_strains.barh(good_species_list_sorted_pretty, proprtion_richness_1, height=0.8, align='center', label='1', color='lightblue')
ax_strains.barh(good_species_list_sorted_pretty, proprtion_richness_2, height=0.8, align='center', label='2', color='dodgerblue', left=proprtion_richness_1)
ax_strains.barh(good_species_list_sorted_pretty, proprtion_richness_3, height=0.8, align='center', label='3', color='royalblue', left=proprtion_richness_1+proprtion_richness_2)
ax_strains.barh(good_species_list_sorted_pretty, proprtion_richness_4,height=0.8, align='center', label='4', color='darkblue', left=proprtion_richness_1+proprtion_richness_2+proprtion_richness_3)


ax_strains.set_xlim([0, 1])
ax_strains.set_ylim([-0.6, len(good_species_list_sorted_pretty)-0.3])

leg = ax_strains.legend(loc="lower left", fontsize=6)
leg.set_title('# strains',prop={'size':7})

ax_strains.set_xlabel('Proportion of hosts', fontsize=11)

ax_strains.xaxis.set_tick_params(labelsize=8)
ax_strains.yaxis.set_tick_params(labelsize=7)

#ax_strains.set_yticklabels(good_species_list_sorted_pretty, rotation=30, ha='right')


fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.15)
fig.savefig("%sstrain_richness_rank_curve.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
