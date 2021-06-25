from __future__ import division
import os, sys
import bz2
import random
import itertools
import config
import parse_midas_data
import numpy
import pickle

import gzip

import matplotlib.pyplot as plt


random.seed(123456789)

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

variant_color_dict = {'1D':'r', '4D':'b'}

good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = ['Bacteroides_vulgatus_57955']

fig, ax = plt.subplots(figsize = (4, 4))

for species_name in good_species_list:

    sys.stderr.write("%s\n" % species_name)

    filename_ = "%sjaccard_distance_snps/%s.pickle" % (data_directory, species_name)

    # Load data (deserialize)
    with open(filename_, 'rb') as handle:
        snp_dict = pickle.load(handle)

    #for allowed_variant_type in allowed_variant_types:

    rhos_1D_list = snp_dict['1D']['jaccard_distances']
    rhos_1D_array = numpy.asarray(rhos_1D_list)

    rhos_4D_list = snp_dict['4D']['jaccard_distances']
    rhos_4D_array = numpy.asarray(rhos_4D_list)

    if (len(rhos_1D_array) < 1000) or (len(rhos_4D_array) < 1000):
        continue

    rho_range = numpy.linspace(0, 1, num=10000)

    rhos_1D_survival = [len(rhos_1D_array[rhos_1D_array >= i]) / float(len(rhos_1D_array)) for i in rho_range]
    rhos_1D_survival = numpy.asarray(rhos_1D_survival)
    #rhos_1D_survival = 1 - rhos_1D_cdf

    rhos_4D_survival = [len(rhos_4D_array[rhos_4D_array >= i]) / float(len(rhos_4D_array)) for i in rho_range]
    rhos_4D_survival = numpy.asarray(rhos_4D_survival)
    #rhos_1D_survival = 1 - rhos_1D_cdf
    #print(rhos_4D_survival)
    difference = rhos_1D_survival - rhos_4D_survival

    #print(rhos_1D_survival)

    #rhos_array_sort = numpy.sort(rhos_array)
    #survival = 1-  numpy.arange(len(rhos_array_sort))/float(len(rhos_array_sort))
    #ax.plot(rho_range, difference, c=variant_color_dict[allowed_variant_type], alpha=0.8)

    ax.plot(rho_range, difference, c='k', alpha=0.6, zorder=2)


ax.axhline(y=0, color='grey', linestyle=':', lw=2, alpha = 1, zorder=1)


#ax.set_xlabel('Correlation beteween pair of SNPs, ' + r'$\rho$' , fontsize = 12)
ax.set_xlabel('Jaccard distance b/w pair of SNPs', fontsize = 12)

#ax.set_ylabel('Difference in prop. SNP pairs ' + r'$\geq \rho$' + '\n Nonsynonymous ' + r'$-$'  + ' synonymous', fontsize = 12)
ax.set_ylabel('Difference in prop. SNP pairs with distance >=  ', fontsize = 12)


fig.subplots_adjust(hspace=0.4, wspace=0.35) #hspace=0.3, wspace=0.5
fig_name =  "%scorrelation_SNP_pairs.png" % parse_midas_data.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
