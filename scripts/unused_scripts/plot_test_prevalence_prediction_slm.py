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
import plot_utils
import stats_utils

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from scipy.stats import gamma, gaussian_kde, ks_2samp

import calculate_predicted_prevalence_mapgd



data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])


clade_type = 'all'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'
#max_n_occurances = 7
species_name = 'Alistipes_finegoldii_56071'


species_color_map, ordered_species_list = plot_utils.get_species_color_map()

#prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()

print(clade_type)

prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict(species_name, clade_type)
species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()


#print(prevalence_dict_mapgd)



#species_name = 'Bacteroides_ovatus_58035'

#print(prevalence_dict_mapgd.keys())

#print(prevalence_dict_mapgd[species_name].keys())

#print(prevalence_dict_mapgd[species_name][clade_type].keys())

#print(prevalence_dict_mapgd[species_name][clade_type][pi_type].keys())



observed_prevalence_mapgd_slm = prevalence_dict_mapgd[pi_type][variant_type]['observed_prevalence_mapgd_slm']
predicted_prevalence_mapgd_slm = prevalence_dict_mapgd[pi_type][variant_type]['predicted_prevalence_mapgd_slm']

#observed_prevalence_mapgd = prevalence_dict_mapgd[pi_type][variant_type]['observed_prevalence_mapgd']
#predicted_prevalence_mapgd = prevalence_dict_mapgd[pi_type][variant_type]['predicted_prevalence_mapgd']


observed_prevalence_mapgd_slm = numpy.asarray(observed_prevalence_mapgd_slm)
predicted_prevalence_mapgd_slm = numpy.asarray(predicted_prevalence_mapgd_slm)

#observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)
#predicted_prevalence_mapgd = numpy.asarray(predicted_prevalence_mapgd)



#error = numpy.absolute(observed_prevalence_mapgd - predicted_prevalence_mapgd) / observed_prevalence_mapgd
#error_slm = numpy.absolute(observed_prevalence_mapgd_slm - predicted_prevalence_mapgd_slm) / observed_prevalence_mapgd_slm

#error_no_nans = error[(~numpy.isnan(error)) & (~numpy.isnan(error_slm)) ]
#error_slm_no_nans = error_slm[(~numpy.isnan(error)) & (~numpy.isnan(error_slm)) ]


#print(prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['sites'])

#print(len(prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd']))

#observed_prevalence_mapgd_slm = numpy.asarray(observed_prevalence_mapgd)
#predicted_prevalence_mapgd_slm = numpy.asarray(predicted_prevalence_mapgd)

observed_prevalence_mapgd_slm_no_nan = observed_prevalence_mapgd_slm[(~numpy.isnan(observed_prevalence_mapgd_slm)) & (~numpy.isnan(predicted_prevalence_mapgd_slm))]
predicted_prevalence_mapgd_slm_no_nan = predicted_prevalence_mapgd_slm[(~numpy.isnan(observed_prevalence_mapgd_slm)) & (~numpy.isnan(predicted_prevalence_mapgd_slm))]


##observed_prevalence_mapgd_no_nan = observed_prevalence_mapgd[(~numpy.isnan(observed_prevalence_mapgd)) & (~numpy.isnan(predicted_prevalence_mapgd))]
#predicted_prevalence_mapgd_no_nan = predicted_prevalence_mapgd[(~numpy.isnan(observed_prevalence_mapgd)) & (~numpy.isnan(predicted_prevalence_mapgd))]


#observed_prevalence = observed_prevalence_mapgd_slm_no_nan[predicted_prevalence_mapgd_slm_no_nan>0]
#predicted_prevalence = predicted_prevalence_mapgd_slm_no_nan[predicted_prevalence_mapgd_slm_no_nan>0]


#print( numpy.absolute(predicted_prevalence - observed_prevalence) / observed_prevalence )


# write code to do KDE on points with prevalence < 0.8

fig, ax = plt.subplots(figsize=(4,4))





all_ = numpy.concatenate([observed_prevalence_mapgd_slm_no_nan, predicted_prevalence_mapgd_slm_no_nan])

print(len(observed_prevalence_mapgd_slm_no_nan))


# Calculate the point density
xy = numpy.vstack([observed_prevalence_mapgd_slm_no_nan, predicted_prevalence_mapgd_slm_no_nan])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = observed_prevalence_mapgd_slm_no_nan[idx], predicted_prevalence_mapgd_slm_no_nan[idx], z[idx]

max_ = max(all_)
min_ = min(all_)

ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
ax.set_xlim([min_, 1.1])
ax.set_ylim([min_, 1.1])


ax.scatter(x, y, c=z, cmap="Blues", s=30, alpha=0.9, edgecolor='', zorder=1)

ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_title(figure_utils.get_pretty_species_name(species_name)  + ' '+ clade_type, fontsize=12, fontweight='bold', color='k' )




fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# , dpi = 600
fig.savefig("%spredicted_observed_prevalence_slm.png" % (config.analysis_directory), format='png', bbox_inches = "tight", dpi = 600, pad_inches = 0.4)
plt.close()
