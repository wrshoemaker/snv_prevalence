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
from collections import Counter

import sys
import bz2
import os
import config
import numpy
from random import sample, shuffle
#import numpy as np

import pickle

import plot_utils
import snps_utils
import sample_utils
import parse_midas_data


data_dir = config.data_directory


#good_species_list = parse_midas_data.parse_good_species_list()


relative_abundance_path = '%sspecies/relative_abundance.txt.bz2' % data_dir


relative_abundance = bz2.BZ2File(relative_abundance_path)
relative_abundance_dict = {}
for line in relative_abundance:
    if 'species_id' in line:
        continue
    line = line.strip().split('\t')
    species = line[0]
    
    relative_abund_i = [float(x) for x in line[1:]]
    relative_abund_i = numpy.asarray(relative_abund_i)

    relative_abundance_dict[species] = {}
    relative_abundance_dict[species]['mean_rel_abund'] = numpy.mean(relative_abund_i)
    relative_abund_i_no_zeros = relative_abund_i[relative_abund_i>0]
    if len(relative_abund_i_no_zeros) > 1:
        relative_abundance_dict[species]['mean_rel_abund_no_zeros'] = numpy.mean(relative_abund_i_no_zeros)

relative_abundance.close()



kegg_dict_all_species = {}

for species in good_species_list:

	kegg_path = '%skegg_pi/kegg_pi_%s.dat' % (parse_midas_data.data_directory, species)

	if os.path.exists(kegg_path) == False:
		continue

	kegg_dict = pickle.load( open(kegg_path, "rb" ) )

	#pathways = kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys()

	for pathway in kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys():

		if (pathway == '') or (pathway == 'Annotated pathways'):
			continue

		pathway_array = kegg_dict['fraction_nonsynonymous_per_pathway_core'][pathway]

		n_hosts = pathway_array.shape[0]

		values = pathway_array[numpy.triu_indices(n_hosts, 1)]

		values = values[(values>0) & (values<1)]

		values = values / (1-values)

		#if pathway not in kegg_dict_all_species:
		#	kegg_dict_all_species[pathway] = []

		if len(values[(values>0) & (values<1)]) <= 20:
			continue
		#kegg_dict_all_species[pathway].append(numpy.mean(values[(values>0) & (values<1)]))

		if pathway not in kegg_dict_all_species:
			kegg_dict_all_species[pathway] = {}
			kegg_dict_all_species[pathway]['mean_dnds_list'] = []
			kegg_dict_all_species[pathway]['species_list'] = []

		kegg_dict_all_species[pathway]['mean_dnds_list'].append(numpy.mean(values))
		kegg_dict_all_species[pathway]['species_list'].append(species)



percentile_dict = {}

for pathway, pathway_dict in kegg_dict_all_species.items():

	mean_dnds_list_pathway = kegg_dict_all_species[pathway]['mean_dnds_list']

	if len(mean_dnds_list_pathway) <10:
		continue

	for species_j, mean_dnds_j in zip(kegg_dict_all_species[pathway]['species_list'], mean_dnds_list_pathway):

		if species_j not in percentile_dict:
			percentile_dict[species_j] = {}
			percentile_dict[species_j]['percentiles_list'] = []
			percentile_dict[species_j]['pathways'] = []

		percentile_dict[species_j]['percentiles_list'].append(mean_dnds_j)

		percentile_dict[species_j]['pathways'].append(pathway)



fig, ax = plt.subplots(figsize=(4,4))

x_list = []
y_list = []

for species in percentile_dict.keys():

    mean_dnds = numpy.mean(percentile_dict[species]['percentiles_list'])

    if species not in relative_abundance_dict:
        continue

    mean_rel_abund = relative_abundance_dict[species]['mean_rel_abund']

    print(species, mean_dnds, mean_rel_abund)

    ax.scatter(mean_rel_abund, mean_dnds)

    x_list.append(mean_rel_abund)
    y_list.append(mean_dnds)


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(x_list), numpy.log10(y_list))

x_log10_fit_range =  numpy.linspace(numpy.log10(min(x_list) * 0.5), numpy.log10(max(x_list) * 1.5), 10000)

y_fit_range = 10 ** (slope*x_log10_fit_range + intercept)

print(p_value)

ax.plot(10**x_log10_fit_range, y_fit_range, c='k', lw=3, linestyle='--', label='OLS fit', zorder=2)


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel('Mean relative abundance', fontsize=10)
ax.set_ylabel('Mean dN/dS', fontsize=10)



fig.subplots_adjust(hspace=0.15, wspace=0.15)
fig.savefig("%s%s.png" % (config.analysis_directory, 'rel_abund_vs_dnds'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
