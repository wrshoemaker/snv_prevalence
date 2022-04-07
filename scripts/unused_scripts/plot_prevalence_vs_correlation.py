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
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_prevalence_mapgd



species_name = 'Bacteroides_xylanisolvens_57185'
clade_type = 'all'
variant_type = '4D'

frequency_dict = calculate_predicted_prevalence_mapgd.load_frequency_dict(species_name, clade_type)
frequency_dict = frequency_dict[variant_type]
# sort sites by number non-zero occurances

sites = list(frequency_dict.keys())
site_pairs = list(combinations(sites, 2))

site_pairs = site_pairs[:10000]

#print(len(sites), len(site_pairs))

r2_dict = {}

proprtion_intersecting_all = []
r2_all = []
for site_pair in site_pairs:

    site_1 = site_pair[0]
    site_2 = site_pair[1]

    site_1_dict = frequency_dict[site_1]
    site_2_dict = frequency_dict[site_2]

    site_1_samples = numpy.asarray(site_1_dict['samples'])
    site_2_samples = numpy.asarray(site_2_dict['samples'])

    site_1_f = numpy.asarray(site_1_dict['frequencies'])
    site_2_f = numpy.asarray(site_2_dict['frequencies'])

    site_intersection = numpy.intersect1d(site_1_samples, site_2_samples)

    if site_intersection < 10:
        continue

    idx_1 = numpy.asarray([numpy.where(site_1_samples == s)[0][0] for s in site_intersection])
    idx_2 = numpy.asarray([numpy.where(site_2_samples == s)[0][0] for s in site_intersection])

    if (len(idx_1) < 10) or (len(idx_2) < 10):
        continue

    site_1_f = site_1_f[idx_1]
    site_2_f = site_2_f[idx_2]

    # keep hosts where at least one site is non zero
    idx_non_zero = (site_1_f>0) | (site_2_f>0)
    if len(idx_non_zero) < 10:
        continue

    site_1_f_to_keep = site_1_f[idx_non_zero]
    site_2_f_to_keep = site_2_f[idx_non_zero]

    n_non_zero_1 = sum(site_1_f_to_keep>0)
    n_non_zero_2 = sum(site_2_f_to_keep>0)

    # only look at pairs with same number of non-zero values
    #if n_non_zero_1 != n_non_zero_2:
    #    continue

    #only keep samples with at least three non-zero values
    if (n_non_zero_1 <= 3):
        continue

    r2 = numpy.corrcoef(site_1_f_to_keep, site_2_f_to_keep)[0,1]

    if numpy.isnan(r2) == True:
        continue

    #print(n_non_zero_1, r2)

    # number of hosts with f > 0 for both snvs
    n_intersecting = sum((site_1_f_to_keep>0) & (site_2_f_to_keep>0))

    proprtion_intersecting = n_intersecting/len(site_1_f_to_keep)


    if n_intersecting not in r2_dict:
        r2_dict[n_intersecting] = []

    #r2_dict[n_intersecting].append(r2)

    #print(proprtion_intersecting, len(site_1_f_to_keep), r2)

    proprtion_intersecting_all.append(proprtion_intersecting)
    r2_all.append(r2)

    #if n_non_zero_1 == n_non_zero_2:

    #    if n_non_zero_1 not in r2_dict:
    #        r2_dict[n_non_zero_1] = []

    #    #print(n_non_zero_1, n_non_zero_2, r2)

    #    r2_dict[n_non_zero_1].append(r2)


#for key, value in r2_dict.items():

#    print(key, len(value), numpy.mean(value))


fig, ax = plt.subplots(figsize=(4,4))

ax.scatter(proprtion_intersecting_all, r2_all, alpha=0.2)


fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
fig.savefig("%scoprevalence_vs_correlation.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
plt.close()
