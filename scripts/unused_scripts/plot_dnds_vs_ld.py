from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
#from mpl_toolkits.axes_grid.inset_locator import inset_axes
import pylab

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from numpy.random import randint, binomial, choice, poisson
from scipy.stats import poisson as poisson_distribution
import scipy.stats as stats

from collections import defaultdict
import sys
import os
import config
import numpy
from random import sample, shuffle
#import numpy as np
import parse_HMP_data
import calculate_linkage_disequilibria
import ld_utils

import pickle

import plot_utils
import snps_utils
import sample_utils
import parse_midas_data

from matplotlib.patches import Rectangle


num_bootstraps = 100
#num_bootstraps = 2
n_percentile_permutations = 10000
#n_percentile_permutations = 10

sub_plot_labels = ['a','b','c']

bin_width_exponent = 1.1
variant_types = ['4D','1D']

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize



data_dir = config.data_directory
picles_dir = "%spickles/" % data_dir
picles_dir = "%spickles/" % data_dir

#syn_differences.pkl

syn_differences = pickle.load(open("%sbetween_syn_differences.pkl" % (picles_dir), 'rb'))
syn_opportunities = pickle.load(open("%sbetween_syn_opportunities.pkl" % (picles_dir), 'rb'))
syn_pseudocounts = pickle.load(open("%sbetween_syn_pseudocounts.pkl" % (picles_dir), 'rb'))
non_differences = pickle.load(open("%sbetween_non_differences.pkl" % (picles_dir), 'rb'))
non_pseudocounts = pickle.load(open("%sbetween_non_pseudocounts.pkl" % (picles_dir), 'rb'))
non_opportunities = pickle.load(open("%sbetween_non_opportunities.pkl" % (picles_dir), 'rb'))


species_color_map, ordered_species_list = plot_utils.get_species_color_map()


def tp_to_category(tp_pair):
		tpa, tpb = tp_pair
		string = 'MI' if (tpa[0], tpb[0]) == ('I', 'M') else tpa[0]+tpb[0]
		return string

# Simplified version

pSs_by_tp_cat = defaultdict(list)
pNpSs_by_tp_cat = defaultdict(list)
pSs_by_species_tp_cat = defaultdict(dict)
pNpSs_by_species_tp_cat = defaultdict(dict)

all_syn_differences_by_tp_cat = defaultdict(list)
all_syn_opportunities_by_tp_cat = defaultdict(list)
all_non_differences_by_tp_cat = defaultdict(list)
all_non_opportunities_by_tp_cat = defaultdict(list)
#all_core_differences_by_tp_cat = defaultdict(list)
all_core_opportunities_by_tp_cat = defaultdict(list)

for tp_pair in syn_differences:
	try:
		tp_cat = tp_to_category(tp_pair)
	except:
		continue

	for species_name in syn_differences[tp_pair]:

		syn_opps = syn_opportunities[tp_pair][species_name]
		syn_diffs = syn_differences[tp_pair][species_name]
		non_opps = non_opportunities[tp_pair][species_name]
		non_diffs = non_differences[tp_pair][species_name]

		pSs = syn_diffs*1.0/syn_opps
		pNs = non_diffs*1.0/non_opps
		pseudo_pSs = 1.0/(syn_opps/2.0+non_opps)
		pseudo_pNs = 1.0/(syn_opps/2.0+non_opps)

		pNpSs = ((pseudo_pNs+pNs)/(pseudo_pSs+pSs))

		good_idxs = ((syn_diffs+non_diffs)>=10)
		pSs_by_tp_cat[tp_cat] += list(pSs[good_idxs])
		pNpSs_by_tp_cat[tp_cat] += list(pNpSs[good_idxs])

		all_syn_differences_by_tp_cat[tp_cat].extend(syn_diffs[good_idxs])
		all_syn_opportunities_by_tp_cat[tp_cat].extend(syn_opps[good_idxs])
		all_non_differences_by_tp_cat[tp_cat].extend(non_diffs[good_idxs])
		all_non_opportunities_by_tp_cat[tp_cat].extend(non_opps[good_idxs])

		try:
			pSs_by_species_tp_cat[species_name][tp_cat] += list(pSs[good_idxs])
			pNpSs_by_species_tp_cat[species_name][tp_cat] += list(pNpSs[good_idxs])
		except:
			pSs_by_species_tp_cat[species_name][tp_cat] = list(pSs[good_idxs])
			pNpSs_by_species_tp_cat[species_name][tp_cat] = list(pNpSs[good_idxs])

median_pSs_by_tp_cat = defaultdict(list)
median_pNpSs_by_tp_cat = defaultdict(list)

pNpSs_species_dict = {}

for species_name in pNpSs_by_species_tp_cat:
    for tp_cat in pNpSs_by_species_tp_cat[species_name]:

        median_pS = numpy.median(pSs_by_species_tp_cat[species_name][tp_cat])
        median_pSs_by_tp_cat[tp_cat].append(median_pS)

        median_pNpS = numpy.median(pNpSs_by_species_tp_cat[species_name][tp_cat])
        mean_pNpS = numpy.mean(pNpSs_by_species_tp_cat[species_name][tp_cat])
        if median_pNpS > 0.2:
            print(tp_cat + ": " + species_name + " has median pNpS " + str(median_pNpS))
        median_pNpSs_by_tp_cat[tp_cat].append(median_pNpS)



        pNpSs_species_dict[species_name] = {}
        pNpSs_species_dict[species_name]['dNdS_median'] = median_pNpS
        pNpSs_species_dict[species_name]['dNdS_mean'] = mean_pNpS




# calculate LD

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species_list = parse_midas_data.parse_good_species_list()

ld_dict_all_species = calculate_linkage_disequilibria.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)

ld_species = list(ld_dict_all_species.keys())

n_cols = 6
species_nested_list  = list(ld_utils.chunks(ld_species, n_cols))
n_rows = len(species_nested_list)

dmin=1
dmax=3000

for species in ld_species:

    ld_dict = ld_dict_all_species[species]

    ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, control_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(ld_dict)

    ld_1D_4D_ratio_filter_rsquareds_dmin_dmax = ld_1D_4D_ratio_filter_rsquareds[(ld_filter_distances>dmin) & (ld_filter_distances<dmax)]

    mean_ld_1D_4D_ratio_filter_rsquareds = numpy.mean(ld_1D_4D_ratio_filter_rsquareds_dmin_dmax)
    median_ld_1D_4D_ratio_filter_rsquareds = numpy.median(ld_1D_4D_ratio_filter_rsquareds_dmin_dmax)

    pNpSs_species_dict[species]['ld_ratio_mean'] = mean_ld_1D_4D_ratio_filter_rsquareds
    pNpSs_species_dict[species]['ld_ratio_median'] = median_ld_1D_4D_ratio_filter_rsquareds



x_list = []
y_list = []

for species in pNpSs_species_dict.keys():

    if ('ld_ratio_mean' not in pNpSs_species_dict[species]) or ('dNdS_mean' not in pNpSs_species_dict[species]):
        continue

    x_list.append(pNpSs_species_dict[species]['dNdS_mean'])
    y_list.append(pNpSs_species_dict[species]['ld_ratio_mean'])



fig, ax = plt.subplots(figsize=(4,4))


ax.scatter(x_list, y_list)



ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel('Mean dN/dS',  fontsize = 10)
ax.set_ylabel('Mean LD ratio',  fontsize = 10)


fig.subplots_adjust(hspace=0.15, wspace=0.15)
fig.savefig("%s%s.png" % (config.analysis_directory, 'dnds_vs_ld_ratio'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
