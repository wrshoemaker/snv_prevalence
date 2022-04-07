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

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_prevalence_mapgd



clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'



test_path = "%spredicted_prevalence_dict_mapgd/Ruminococcus_bicirculans_59300_all.dat" % (config.data_directory)

with open(test_path, 'rb') as handle:
    prevalence_dict_mapgd = pickle.load(handle)
    #return b



fig, ax = plt.subplots(figsize=(4,4))

predicted_prevalence = prevalence_dict_mapgd[pi_type][variant_type]['predicted_prevalence_mapgd_slm']
predicted_prevalence = numpy.asarray(predicted_prevalence)

observed_prevalence = prevalence_dict_mapgd[pi_type][variant_type]['observed_prevalence_mapgd_slm']
observed_prevalence = numpy.asarray(observed_prevalence)

#print(prevalence_dict_mapgd[pi_type][variant_type].keys())


predicted_prevalence_mapgd = prevalence_dict_mapgd[pi_type][variant_type]['predicted_prevalence_mapgd_slm_best']
predicted_prevalence_mapgd = numpy.asarray(predicted_prevalence_mapgd)

observed_prevalence_mapgd = prevalence_dict_mapgd[pi_type][variant_type]['observed_prevalence_mapgd_slm_best']
observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)





#predicted_prevalence_no_zeros = predicted_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]
#observed_prevalence_no_zeros = observed_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]


predicted_prevalence_no_zeros = predicted_prevalence_mapgd[(observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0) ]
observed_prevalence_no_zeros = observed_prevalence_mapgd[(observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0) ]


#mre = numpy.mean(numpy.absolute(predicted_prevalence_mapgd_no_zeros-observed_prevalence_mapgd_no_zeros)/observed_prevalence_mapgd_no_zeros)

#mre_best = numpy.mean(numpy.absolute(predicted_prevalence_no_zeros-observed_prevalence_no_zeros)/observed_prevalence_no_zeros)

#print("adjusted estimate", mre)
#print("original estimate", mre_best)


all_ = numpy.concatenate([predicted_prevalence_no_zeros,observed_prevalence_no_zeros])

xy = numpy.vstack([observed_prevalence_no_zeros, predicted_prevalence_no_zeros])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = observed_prevalence_no_zeros[idx], predicted_prevalence_no_zeros[idx], z[idx]
ax.scatter(x, y, c=z, cmap="Blues", s=90, alpha=0.9, edgecolor='', zorder=1)


#ax.scatter(observed_prevalence, predicted_prevalence, c='dodgerblue', s=90, alpha=0.9, edgecolor='', zorder=1)

max_ = max(all_)*1.1
min_ = min(all_)*0.8

#print(min(f_max_no_zeros))

ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
ax.set_xlim([min_, max_])
ax.set_ylim([min_, max_])


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)



fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# dpi = 600
fig.savefig("%spredicted_observed_prevalence_mapgd_test.png" % config.analysis_directory, format='png', bbox_inches = "tight", dpi=600, pad_inches = 0.4)
plt.close()
