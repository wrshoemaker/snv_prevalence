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
import parse_midas_data
import prevalence_utils
import scipy.stats as stats

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

import calculate_predicted_prevalence_mapgd
import calculate_predicted_prevalence


prevalence_dict = calculate_predicted_prevalence.load_predicted_prevalence_subsample_dict()
prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict_mapgd.keys())

clade_type = 'all'





fig, ax = plt.subplots(figsize=(4,4))

mre_all = []
mre_mapgd_all = []

for species_name in species_list:

    print(prevalence_dict_mapgd[species_name]['4D'].keys())

    predicted_prevalence_mapgd = prevalence_dict_mapgd[species_name]['4D']['predicted_prevalence']
    predicted_prevalence_mapgd = numpy.asarray(predicted_prevalence_mapgd)

    observed_prevalence_mapgd = prevalence_dict_mapgd[species_name]['4D']['observed_prevalence']
    observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)

    if len(predicted_prevalence_mapgd) < 30:
        continue

    predicted_prevalence_mapgd_no_zeros = predicted_prevalence_mapgd[(observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0)]
    observed_prevalence_mapgd_no_zeros = observed_prevalence_mapgd[(observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0)]

    mre_mapgd = numpy.mean(numpy.absolute(observed_prevalence_mapgd_no_zeros - predicted_prevalence_mapgd_no_zeros) / observed_prevalence_mapgd_no_zeros)


    mre = prevalence_dict[species_name][clade_type]['4D']['MRE']

    predicted_prevalence = prevalence_dict[species_name][clade_type]['4D']['predicted_prevalence_to_plot']
    predicted_prevalence = numpy.asarray(predicted_prevalence)

    observed_prevalence = prevalence_dict[species_name][clade_type]['4D']['observed_prevalence_to_plot']
    observed_prevalence = numpy.asarray(observed_prevalence)

    predicted_prevalence_subset = predicted_prevalence[observed_prevalence <= max(observed_prevalence_mapgd_no_zeros)]
    observed_prevalence_subset = observed_prevalence[observed_prevalence <= max(observed_prevalence_mapgd_no_zeros)]

    mre_subset = numpy.mean(numpy.absolute(observed_prevalence_subset - predicted_prevalence_subset) / observed_prevalence_subset)


    mre_all.append(mre)
    mre_mapgd_all.append(mre_mapgd)


ax.scatter(mre_all, mre_mapgd_all, c='k', alpha=0.7, label='Species', zorder=2)



slope, intercept, r_value, p_value, std_err = stats.linregress(mre_all, numpy.log10(mre_mapgd_all))
x_range, y_pred, lcb, ucb = prevalence_utils.get_confidence_hull(mre_all, numpy.log10(mre_mapgd_all))

print(p_value)

#ax.plot(x_range, 10**y_pred, color='k', linestyle='--', linewidth=2, zorder=3, label='OLS regression')
#ax.plot(x_range, 10**lcb, color='k', linestyle=':', linewidth=2, zorder=3, label=r'$95\%$' + ' Confidence band')
#ax.plot(x_range, 10**ucb, color='k', linestyle=':', linewidth=2, zorder=3)

#ax.legend(loc="lower right", prop={'size': 6})

min_ = min(mre_all)*0.8
max_ = max(mre_all)*1.1

ax.set_xlim([min(mre_all)*0.8, max(mre_all)*1.1])
ax.set_ylim([min(mre_mapgd_all)*0.8, max(mre_mapgd_all)*1.1])

#ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel('Mean relative error, naive', fontsize=12)
ax.set_ylabel('Mean relative error, MAPGD', fontsize=12)

#ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=1)


fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%serror_naive_vs_mapgd.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
