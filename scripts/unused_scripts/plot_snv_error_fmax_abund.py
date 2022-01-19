from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import scipy.stats as stats


import diversity_utils
import parse_midas_data
import calculate_predicted_occupancy
import prevalence_utils

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import figure_utils


data_dir = config.data_directory


data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict.keys())
species_list.sort()

good_species = 'Eubacterium_rectale_56927'
#bad_species = 'Bacteroides_vulgatus_57955'



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


    if species_name in [good_species]:
        good_bad_species_value_dict[species_name] = {}
        good_bad_species_value_dict[species_name]['mre'] = mre
        good_bad_species_value_dict[species_name]['mean_relative_abundance'] = mean_relative_abund_i

relative_abundance.close()

#print(error_array)

error_array = numpy.asarray(error_array)

print(sum(error_array<0.5), len(error_array))

print(max(error_array))




fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.15)

ax_f_max = plt.subplot2grid((1, 2), (0,0))
ax_abundance = plt.subplot2grid((1, 2), (0,1))


fontsize=8
ax_f_max.text(0.01, 1.07, figure_utils.sub_plot_labels[0], fontsize=fontsize, fontweight='bold', ha='center', va='center', transform=ax_f_max.transAxes)
ax_abundance.text(0.01, 1.07, figure_utils.sub_plot_labels[1], fontsize=fontsize, fontweight='bold', ha='center', va='center', transform=ax_abundance.transAxes)




good_bad_color_dict = prevalence_utils.good_bad_color_dict

for species_name in species_list:

    f_max_line = prevalence_dict[species_name]['all']['4D']['f_max_line']
    f_max_line = numpy.asarray(f_max_line)


    f_max_relative_error_line = prevalence_dict[species_name]['all']['4D']['f_max_vs_relative_error_line']


    if species_name == good_species:
        c = good_bad_color_dict[species_name]
        alpha=1
        zorder=2
        lw=3

        ax_f_max.plot(10**f_max_line, f_max_relative_error_line, ls ='-', c=c, lw=lw, alpha=alpha, label=figure_utils.get_pretty_species_name(species_name), zorder=zorder)

    else:
        c = 'k'
        alpha=0.4
        zorder=1
        lw=2

        ax_f_max.plot(10**f_max_line, f_max_relative_error_line, ls ='-', c=c, lw=lw, alpha=alpha, zorder=zorder)





ax_f_max.axvline(x=0.2, ls='--', lw=2, c='k')
ax_f_max.axvline(x=0.8, ls='--', lw=2, c='k')

ax_f_max.set_xscale('log', basex=10)
ax_f_max.set_xlabel('Maximum frequency\nacross hosts, ' + r'$f_{max}$', fontsize=14)
ax_f_max.set_ylabel('SNV prevalence MRE for an ' + r'$f_{max}$' + ' bin', fontsize=12)



ax_f_max.legend(loc="upper left", prop={'size': 6})


ax_abundance.scatter(relative_abundance_array, error_array, c='k', alpha=0.7)


# plot good and bad species
for species_name in [good_species]:

    ax_abundance.scatter(good_bad_species_value_dict[species_name]['mean_relative_abundance'], good_bad_species_value_dict[species_name]['mre'], c=good_bad_color_dict[species_name], alpha=1, s=60, label=figure_utils.get_pretty_species_name(species_name), zorder=4)



ax_abundance.set_xlim([min(relative_abundance_array)*0.8, max(relative_abundance_array)*1.1])
ax_abundance.set_ylim([min(error_array)*0.8, max(error_array)*1.1])


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(relative_abundance_array), error_array)


x_range, y_pred, lcb, ucb = prevalence_utils.get_confidence_hull(numpy.log10(relative_abundance_array), error_array)


ax_abundance.text(0.5,0.95, r'$r^{2} = $' + str(round(r_value**2, 3)), fontsize=8, color='k', ha='center', va='center', transform=ax_abundance.transAxes  )
ax_abundance.text(0.522,0.885, r'$P = $' + str(round(p_value, 5)), fontsize=8, color='k', ha='center', va='center', transform=ax_abundance.transAxes  )



ax_abundance.plot(10**x_range, y_pred, color='k', linestyle='--', linewidth=2, zorder=3, label='OLS regression')

ax_abundance.plot(10**x_range, lcb, color='k', linestyle=':', linewidth=2, zorder=3, label=r'$95\%$' + ' Confidence band')
ax_abundance.plot(10**x_range, ucb, color='k', linestyle=':', linewidth=2, zorder=3)


ax_abundance.legend(loc="upper left", prop={'size': 6})



ax_abundance.set_xscale('log', basex=10)
#ax.set_yscale('log', basey=10)


ax_abundance.set_xlabel('Mean relative species abundance', fontsize=12)
ax_abundance.set_ylabel('SNV prevalence MRE', fontsize=12)

#ax.set_ylim([0.1, 1])

fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%serror_fmax_abund.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
