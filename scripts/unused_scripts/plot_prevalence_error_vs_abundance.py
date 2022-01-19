from __future__ import division
import sample_utils
import config
import parse_midas_data
import os.path
import os
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
from math import log10,ceil,fabs,isnan
from numpy.random import randint, choice, multinomial

import prevalence_utils

import midas_db_utils

import scipy.stats as stats

import parse_HMP_data
import figure_utils
import calculate_predicted_occupancy

#good_species_list = [good_species_list[3]]

import matplotlib.pyplot as plt
import matplotlib.cm as cm

data_dir = config.data_directory


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()


relative_abundance_path = '%sspecies/relative_abundance.txt.bz2' % data_dir


relative_abundance_array = []
error_array = []
relative_abundance = bz2.BZ2File(relative_abundance_path)
#relative_abundance_dict = {}
samples = relative_abundance.readline()
samples = samples.strip().split('\t')
samples = [str(x) for x in samples[1:]]

good_bad_species_value_dict = {}

for line_idx, line in enumerate(relative_abundance):
    line = line.strip().split('\t')
    species_name = line[0]

    if species_name not in prevalence_dict:
        continue

    samples_species = prevalence_dict[species_name]['all']['samples']

    # get idxs of samples

    samples_species_idx =  [samples.index(s) for s in samples_species]
    samples_species_idx = numpy.asarray(samples_species_idx)


    relative_abund_i = [float(x) for x in line[1:]]
    relative_abund_i = numpy.asarray(relative_abund_i)

    # add code to only include samples used in the analysis

    mean_relative_abund_i = numpy.mean(relative_abund_i[samples_species_idx])

    relative_abundance_array.append(mean_relative_abund_i)

    mre = prevalence_dict[species_name]['all']['4D']['MRE']

    error_array.append(mre)


    if species_name in prevalence_utils.good_bad_color_dict:
        good_bad_species_value_dict[species_name] = {}
        good_bad_species_value_dict[species_name]['mre'] = mre
        good_bad_species_value_dict[species_name]['mean_relative_abundance'] = mean_relative_abund_i

relative_abundance.close()




fig, ax = plt.subplots(figsize=(4,4))


ax.scatter(relative_abundance_array, error_array, c='k', alpha=0.7)


# plot good and bad species
for species_name in prevalence_utils.good_bad_color_dict.keys():

    ax.scatter(good_bad_species_value_dict[species_name]['mean_relative_abundance'], good_bad_species_value_dict[species_name]['mre'], c=prevalence_utils.good_bad_color_dict[species_name], alpha=1, s=60, label=figure_utils.get_pretty_species_name(species_name), zorder=4)





ax.set_xlim([min(relative_abundance_array)*0.8, max(relative_abundance_array)*1.1])
ax.set_ylim([min(error_array)*0.8, max(error_array)*1.1])


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(relative_abundance_array), error_array)


x_range, y_pred, lcb, ucb = prevalence_utils.get_confidence_hull(numpy.log10(relative_abundance_array), error_array)


ax.text(0.5,0.95, r'$r^{2} = $' + str(round(r_value**2, 3)), fontsize=8, color='k', ha='center', va='center', transform=ax.transAxes  )
ax.text(0.522,0.885, r'$P = $' + str(round(p_value, 5)), fontsize=8, color='k', ha='center', va='center', transform=ax.transAxes  )



ax.plot(10**x_range, y_pred, color='k', linestyle='--', linewidth=2, zorder=3, label='OLS regression')

ax.plot(10**x_range, lcb, color='k', linestyle=':', linewidth=2, zorder=3, label=r'$95\%$' + ' Confidence band')
ax.plot(10**x_range, ucb, color='k', linestyle=':', linewidth=2, zorder=3)


ax.legend(loc="upper left", prop={'size': 6})



ax.set_xscale('log', basex=10)
#ax.set_yscale('log', basey=10)


ax.set_xlabel('Mean relative abundance', fontsize=12)
ax.set_ylabel('Mean relative error', fontsize=12)

#ax.set_ylim([0.1, 1])

fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%smean_abundance_vs_error.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
