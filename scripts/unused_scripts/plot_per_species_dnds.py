from __future__ import division
from math import log10,ceil
from numpy.random import randint


from numpy.random import randint, binomial, choice, poisson
import scipy.stats as stats

from collections import defaultdict
import sys
import os
import config
import numpy
from random import sample, shuffle
#import numpy as np

import pickle

#import plot_utils
import snps_utils
import sample_utils
import parse_midas_data



#unzip negative_selection_data.zip "negative_selection_data/linkage_disequilibria/" -d ./negative_selection_data_pickles



num_bootstraps = 100
#num_bootstraps = 2
n_percentile_permutations = 10000
#n_percentile_permutations = 10

sub_plot_labels = ['a','b','c']


#data_dir = config.data_directory
#data_dir = os.path.expanduser("~/GitHub/negative_selection_microbiome/data/")
#picles_dir = "%spickles/" % data_dir
#picles_dir = "%spickles/" % data_dir

picles_dir = '/u/home/w/wrshoema/project-ngarud/negative_selection_data_pickles/'
#syn_differences.pkl

syn_differences = pickle.load(open("%sbetween_syn_differences.pkl" % (picles_dir), 'rb'))
syn_opportunities = pickle.load(open("%sbetween_syn_opportunities.pkl" % (picles_dir), 'rb'))
syn_pseudocounts = pickle.load(open("%sbetween_syn_pseudocounts.pkl" % (picles_dir), 'rb'))
non_differences = pickle.load(open("%sbetween_non_differences.pkl" % (picles_dir), 'rb'))
non_pseudocounts = pickle.load(open("%sbetween_non_pseudocounts.pkl" % (picles_dir), 'rb'))
non_opportunities = pickle.load(open("%sbetween_non_opportunities.pkl" % (picles_dir), 'rb'))


#species_color_map, ordered_species_list = plot_utils.get_species_color_map()


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

for species_name in pNpSs_by_species_tp_cat:
	for tp_cat in pNpSs_by_species_tp_cat[species_name]:

		median_pS = numpy.median(pSs_by_species_tp_cat[species_name][tp_cat])
		median_pSs_by_tp_cat[tp_cat].append(median_pS)

		median_pNpS = numpy.median(pNpSs_by_species_tp_cat[species_name][tp_cat])
		if median_pNpS > 0.2:
			print(tp_cat + ": " + species_name + " has median pNpS " + str(median_pNpS))
		median_pNpSs_by_tp_cat[tp_cat].append(median_pNpS)




# Bootstrapping dN/dS

avg_cf_ratios = defaultdict(list)
std_cf_ratios = defaultdict(list)
median_cf_ratios = defaultdict(list)
lower_cf_ratios = defaultdict(list)
upper_cf_ratios = defaultdict(list)
avg_sf_ratios = defaultdict(list)
std_sf_ratios = defaultdict(list)

ds = numpy.logspace(-5,-2,50)

#for tp_cat in ['II', 'AA', 'MI', 'MM']:
for tp_cat in ['AA']:

	cf_ratios = [] # cumulative estimates <= total d
	sf_ratios = [] # cumulative estimates >= total d

	sys.stderr.write("Bootstrapping dN/dS...\n")

	for bootstrap_idx in range(num_bootstraps):

		lower_pNpSs, upper_pNpSs = [], []

		all_non_diffs = numpy.array(all_non_differences_by_tp_cat[tp_cat])
		all_non_opps = numpy.array(all_non_opportunities_by_tp_cat[tp_cat])
		all_syn_diffs = numpy.array(all_syn_differences_by_tp_cat[tp_cat])
		all_syn_opps = numpy.array(all_syn_opportunities_by_tp_cat[tp_cat])
		all_NS_differences = all_syn_diffs + all_non_diffs

		# Bootstrap dataset using poisson resampling
		# Each observed difference is considered lambda for Poisson
		# distribution. Resample according to pmf in which output N
		# is weighted by probability that that number of differences
		# occurs within interval lambda, the average number of diffs.

		# When lambda is small, Poisson distribution is highly
		# right skewed. As lambda approaches infinity, become
		# more and more like the binomial distribution

		# Pseudocounts so things w/ 0 counts are not "stuck" in resampling
		# Pseudocounts are chosen w/ dN/dS=1, so should be conservative?
		# (alternatively, we could choose dN/dS=0.1 -- unfair?)

		pseudocount = 0 # 1.0
		bs_non_differences = poisson(all_non_diffs + pseudocount)
		bs_syn_differences = poisson(all_syn_diffs + (all_syn_opps*pseudocount/all_non_opps))

		bs_NS_differences = bs_non_differences + bs_syn_differences

		# Cut down numbers by half on average
		#if bs_syn_differences <= 0:
		#	continue
		print(bs_syn_differences)
		bs_thinned_syn_differences_1 = binomial(bs_syn_differences, 0.5)
		bs_thinned_syn_differences_2 = bs_syn_differences - bs_thinned_syn_differences_1

		# Bootstrapped dS
		bs_divergence = bs_thinned_syn_differences_1 / (all_syn_opps/2.0)

		for d in ds:

			lower_idxs = (bs_divergence <= d)*(all_NS_differences>0.5)*(bs_NS_differences>0.5)
			upper_idxs = (bs_divergence > d)*(all_NS_differences>0.5)*(bs_NS_differences>0.5)

			if lower_idxs.sum()<1.5:
					lower_pNpSs.append(-1)
			else:
					lower_cumulative_non_differences = (bs_non_differences)[lower_idxs].sum()
					lower_cumulative_expected_non_differences = (bs_thinned_syn_differences_2[lower_idxs]*2.0/all_syn_opps[lower_idxs]*all_non_opps[lower_idxs]).sum()
					lower_pNpSs.append( (lower_cumulative_non_differences)/(lower_cumulative_expected_non_differences) )

			if upper_idxs.sum()<1.5:
					upper_pNpSs.append(-1)
			else:
					upper_cumulative_non_differences = (bs_non_differences[upper_idxs]).sum()
					upper_cumulative_expected_non_differences = (bs_thinned_syn_differences_2[upper_idxs]*2.0/all_syn_opps[upper_idxs]*all_non_opps[upper_idxs]).sum()
					upper_pNpSs.append( (upper_cumulative_non_differences)/(upper_cumulative_expected_non_differences) )

		cf_ratios.append(lower_pNpSs)
		sf_ratios.append(upper_pNpSs)

		if bootstrap_idx % 10 == 0:
			print("On bootstrap %i..." % bootstrap_idx)

	cf_ratios = numpy.array(cf_ratios)
	sf_ratios = numpy.array(sf_ratios)

	for i in range(len(ds)):

		ratios = numpy.sort(cf_ratios[:,i])
		good_idxs = (ratios>-0.5)
		if good_idxs.sum()<1.5:
				avg_cf_ratios[tp_cat].append(-1)
				std_cf_ratios[tp_cat].append(0)

		else:
				median_cf_ratios[tp_cat].append(numpy.median(ratios[good_idxs]))
				idx = long(0.025*good_idxs.sum())
				lower_cf_ratios[tp_cat].append( ratios[good_idxs][idx] )
				upper_cf_ratios[tp_cat].append(ratios[good_idxs][-idx-1])

				avg_cf_ratios[tp_cat].append( ratios[good_idxs].mean() )
				std_cf_ratios[tp_cat].append( ratios[good_idxs].std() )

		ratios = sf_ratios[:,i]
		good_idxs = (ratios>-0.5)
		if good_idxs.sum()<1.5:
				avg_sf_ratios[tp_cat].append(-1)
				std_sf_ratios[tp_cat].append(0)
		else:
				avg_sf_ratios[tp_cat].append( ratios[good_idxs].mean() )
				std_sf_ratios[tp_cat].append( ratios[good_idxs].std() )

	avg_cf_ratios[tp_cat] = numpy.array(avg_cf_ratios[tp_cat])
	std_cf_ratios[tp_cat] = numpy.array(std_cf_ratios[tp_cat])
	median_cf_ratios[tp_cat] = numpy.array(median_cf_ratios[tp_cat])
	upper_cf_ratios[tp_cat] = numpy.array(upper_cf_ratios[tp_cat])
	lower_cf_ratios[tp_cat] = numpy.array(lower_cf_ratios[tp_cat])
	avg_sf_ratios[tp_cat] = numpy.array(avg_sf_ratios[tp_cat])
	std_sf_ratios[tp_cat] = numpy.array(std_sf_ratios[tp_cat])

#tp_cat_descrip = {'II': 'infant-infant', 'AA': 'adult-adult', 'MI': 'mother-infant', 'MM': 'mother-mother'}
tp_cat_descrip = {'AA': 'adult-adult'}


# Create a species color legend
all_species = set()
for species in pSs_by_species_tp_cat:
	all_species.add(species)

all_species_ordered = []
#colors_ordered = []
#for species in ordered_species_list:
#	if species in all_species:
#		all_species_ordered.append(species)
#		colors_ordered.append(species_color_map[species])



tp_cat = 'AA'


# Moving on...
all_pSs = pSs_by_tp_cat[tp_cat]
all_pNpSs = pNpSs_by_tp_cat[tp_cat]


output_file = open('/u/home/w/wrshoema/project-ngarud/mean_dnds.txt', "w")

for species in pSs_by_species_tp_cat:
	if tp_cat in pSs_by_species_tp_cat[species]:
		pSs = pSs_by_species_tp_cat[species][tp_cat]
		pNpSs = pNpSs_by_species_tp_cat[species][tp_cat]
		#color = species_color_map[species]
		#divergence_axis.loglog(pSs, pNpSs, '.', color=color, markersize=6,alpha=0.6,markeredgewidth=0,zorder=0,rasterized=True)

		# mean
		mean_pNpSs = numpy.mean(pNpSs)

		output_file.write("\t".join([species, str(mean_pNpSs)]))
		output_file.write("\n")

output_file.close()
