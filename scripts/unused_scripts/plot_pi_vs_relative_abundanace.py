from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import diversity_utils
import figure_utils
import plot_utils
import parse_midas_data
import parse_HMP_data

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from scipy.stats import gamma, gaussian_kde
import scipy.stats as stats

import prevalence_utils

#import calculate_predicted_occupancy
import calculate_predicted_prevalence_mapgd

data_dir = config.data_directory


clade_type = 'no_strains'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'
max_n_occurances = 7


prevalence_dict = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()

species_color_map, ordered_species_list = plot_utils.get_species_color_map()

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

    #samples_species = prevalence_dict[species_name][clade_type][pi_type][variant_type]
    #print(samples_species.keys())

    pi_dict = calculate_predicted_prevalence_mapgd.load_pi_dict(species_name)

    print(pi_dict.keys())

    pi_samples = list(pi_dict[variant_type].keys())
    # get idxs of samples

    pi_list = [pi_dict[variant_type][s][pi_type] for s in pi_samples]
    pi_list = numpy.asarray(pi_list)
    print(pi_list)
    #mean_log_pi = numpy.mean(numpy.log10(pi_list))
    mean_log_pi = numpy.mean(pi_list)


    #samples_species_idx =  [samples.index(s) for s in samples_species]
    samples_species_idx =  [samples.index(s) for s in pi_samples]

    samples_species_idx = numpy.asarray(samples_species_idx)


    relative_abund_i = [float(x) for x in line[1:]]
    relative_abund_i = numpy.asarray(relative_abund_i)

    # add code to only include samples used in the analysis

    mean_relative_abund_i = numpy.mean(relative_abund_i[samples_species_idx])

    relative_abundance_array.append(mean_relative_abund_i)

    #mre = prevalence_dict[species_name]['all']['4D']['MRE']


    error_array.append(mean_log_pi)


    if species_name in prevalence_utils.good_bad_color_dict:
        good_bad_species_value_dict[species_name] = {}
        good_bad_species_value_dict[species_name]['mre'] = mre
        good_bad_species_value_dict[species_name]['mean_relative_abundance'] = mean_relative_abund_i
        good_bad_species_value_dict[species_name]['mean_log10_pi'] = mean_log_pi


relative_abundance.close()



fig, ax = plt.subplots(figsize=(4,4))


ax.scatter(relative_abundance_array, error_array, c='k', alpha=0.7)


# plot good and bad species
for species_name in prevalence_utils.good_bad_color_dict.keys():

    ax.scatter(good_bad_species_value_dict[species_name]['mean_relative_abundance'], good_bad_species_value_dict[species_name]['mean_log10_pi'], c=species_color_map[species_name], alpha=1, s=60, label=figure_utils.get_pretty_species_name(species_name), zorder=4)





ax.set_xlim([min(relative_abundance_array)*0.8, max(relative_abundance_array)*1.1])
ax.set_ylim([min(error_array)*0.8, max(error_array)*1.1])


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(relative_abundance_array), numpy.log10(error_array))
x_range, y_pred, lcb, ucb = prevalence_utils.get_confidence_hull(numpy.log10(relative_abundance_array), error_array)

#ax.text(0.5,0.95, r'$r^{2} = $' + str(round(r_value**2, 3)), fontsize=8, color='k', ha='center', va='center', transform=ax.transAxes  )
#ax.text(0.522,0.885, r'$P = $' + str(round(p_value, 5)), fontsize=8, color='k', ha='center', va='center', transform=ax.transAxes  )



#ax.plot(10**x_range, y_pred, color='k', linestyle='--', linewidth=2, zorder=3, label='OLS regression')
#ax.plot(10**x_range, lcb, color='k', linestyle=':', linewidth=2, zorder=3, label=r'$95\%$' + ' Confidence band')
#ax.plot(10**x_range, ucb, color='k', linestyle=':', linewidth=2, zorder=3)
#ax.legend(loc="upper left", prop={'size': 6})



ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)


ax.set_xlabel('Mean relative abundance', fontsize=12)
ax.set_ylabel('Mean nucleotide diversity, ' + r'$\hat{\theta}$', fontsize=12)

#ax.set_ylim([0.1, 1])

fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%smean_abundance_vs_pi.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
