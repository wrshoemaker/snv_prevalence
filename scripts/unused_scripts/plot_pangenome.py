import  matplotlib.pyplot as plt

import parse_midas_data
#import pylab
#from pylab import *
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
import os
#import pandas
import parse_patric
import pickle
import sample_utils
import parse_HMP_data

import itertools

import scipy.stats as stats



good_species_list = parse_midas_data.parse_good_species_list()
#good_species_list = ["Eubacterium_rectale_56927"]

#kegg_dict_all_species = {}

prevalence_bins = numpy.linspace(0, 1, num=11)
prevalence_bins_pairs = [(prevalence_bins[i], prevalence_bins[i+1]) for i in range(len(prevalence_bins)-1)]
prevalence_bins_mins = prevalence_bins[:-1]

fig, ax = plt.subplots()


for species_name in good_species_list:

    _path = '%sprevalence_pi/prevalence_pi_%s.dat' % (parse_midas_data.data_directory, species_name)

    if os.path.exists(_path) == False:
        continue

    _dict = pickle.load( open(_path, "rb" ) )

    pi_list = [_dict[x]['avg_pi_matrix_avg'] for x in prevalence_bins_pairs]

    ax.plot(prevalence_bins_mins, pi_list, marker='o')

ax.set_yscale('log', basey=10)



ax.set_xticks(prevalence_bins_mins)
xticklabels = ['%.1f-%.1f' % (prevalence_bins[i], prevalence_bins[i+1]) for i in range(len(prevalence_bins)-1)]

ax.set_xticklabels(xticklabels, fontsize=8)

ax.set_xlabel('Prevalence bin', fontsize=12)
ax.set_ylabel('Average synonymous pairwise divergence', fontsize=12)



fig.subplots_adjust(hspace=0.4, wspace=0.35) #hspace=0.3, wspace=0.5
fig_name =  "%sprevalence_pi.png" % parse_midas_data.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
