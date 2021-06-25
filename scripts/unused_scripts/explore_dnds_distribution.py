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


good_species_list = parse_midas_data.parse_good_species_list()
#good_species_list = ["Eubacterium_rectale_56927"]

kegg_dict_all_species = {}



fig, ax = plt.subplots(figsize=(4,4))


all_pathways = []
for species in good_species_list:

    kegg_path = '%skegg_pi/kegg_pi_%s.dat' % (parse_midas_data.data_directory, species)

    if os.path.exists(kegg_path) == False:
		continue

    kegg_dict = pickle.load( open(kegg_path, "rb" ))

    pathways = kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys()
    all_pathways.extend(pathways)





all_pathways_count = dict(Counter(all_pathways))

pathways_to_keep = [x for x in all_pathways_count.keys() if all_pathways_count[x] > 50 if (x != '') and (x != 'Annotated pathways')]

all_values = []

kegg_dict_all_species = {}

for species in good_species_list:

    kegg_path = '%skegg_pi/kegg_pi_%s.dat' % (parse_midas_data.data_directory, species)

    if os.path.exists(kegg_path) == False:
		continue

    kegg_dict = pickle.load( open(kegg_path, "rb" ))

    #print(kegg_dict)
	#pathways = kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys()
    pathway_means = []
    #all_values = []

    #kegg_dict['dn_ds_per_pathway'] = {}

    for pathway in kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys():

        if (pathway == '') or (pathway == 'Annotated pathways'):
			continue

        if pathway not in pathways_to_keep:
            continue

        if pathway not in kegg_dict_all_species:
            kegg_dict_all_species[pathway] = {}

        pathway_array = kegg_dict['fraction_nonsynonymous_per_pathway_core'][pathway]

        n_hosts = pathway_array.shape[0]

        values = pathway_array[numpy.triu_indices(n_hosts, 1)]

        values = values[(values>0) & (values<1)]

        values = values / (1-values)

        if len(values[(values>0) & (values<1)]) <= 20:
			continue

        if len(values) < 1300:
            continue

        all_values.append(values)

        kegg_dict_all_species[pathway][species] = values

        #kegg_dict['dn_ds_per_pathway'][pathway] = values



all_values = numpy.concatenate(all_values).ravel()
all_values = numpy.log10(all_values)

all_values_std = numpy.std(all_values)
all_values_mean = numpy.mean(all_values)


pathways_to_keep_subset = pathways_to_keep[:5]
fig = plt.figure(figsize = (4*len(pathways_to_keep_subset), 4)) #
fig.subplots_adjust(bottom= 0.15)

for pathway_idx, pathway in enumerate(pathways_to_keep_subset):

    ax = plt.subplot2grid((1, 1*len(pathways_to_keep_subset)), (0, pathway_idx), colspan=1)
    ax.set_title(pathway, fontsize=9, fontweight='bold' )

    mean_all = []
    var_all = []

    for species, values in kegg_dict_all_species[pathway].items():

        values = numpy.log10(values)

        values = (values-all_values_mean) / all_values_std

        counts, bin_edges = numpy.histogram(values, 25, density=True)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.

        ax.scatter(bin_centres, counts)



        #values_mean = numpy.mean(values)
        #values_var = numpy.var(values)

        #mean_all.append(values_mean)
        #var_all.append(values_var)

        #ax.scatter(values_mean, values_var)

    ax.set_yscale('log', basey=10)

    ax.set_ylim([0.001, 3])

    #print(zip(numpy.log10(mean_all), numpy.log10(var_all)))

    #slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_all), numpy.log10(var_all))

    #x_log10_fit_range =  numpy.linspace(numpy.log10(min(mean_all) * 0.5), numpy.log10(max(mean_all) * 1.5), 10000)

    #y_fit_range = 10 ** (slope*x_log10_fit_range + intercept)

    #ax.plot(10**x_log10_fit_range, y_fit_range, c='k', lw=3, linestyle='--', label='OLS fit', zorder=2)

    #ax.set_xscale('log', basex=10)
    #ax.set_yscale('log', basey=10)


    ax.set_xlabel('Standardized log10 dN/dS', fontsize=10)
    ax.set_ylabel('Probability density', fontsize=10)



#ax.set_xscale('log', basex=10)

fig.subplots_adjust(hspace=0.4, wspace=0.4)
fig.savefig("%s%s.png" % (config.analysis_directory, 'dnds_distribution'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
