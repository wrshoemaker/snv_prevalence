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
from matplotlib.lines import Line2D

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_occupancy

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])
clade_types = ['all','largest_clade']


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()


species_list = list(prevalence_dict.keys())

mre_4d = [prevalence_dict[species_name]['all']['4D']['MRE'] for species_name in species_list]
mre_1d = [prevalence_dict[species_name]['all']['1D']['MRE'] for species_name in species_list]



fig, ax = plt.subplots(figsize=(4,4))

ax.plot([min(mre_4d)*0.8, max(mre_4d)*1.1],[min(mre_1d)*0.8, max(mre_1d)*1.1], ls='--', c='k')

ax.scatter(mre_4d, mre_1d, c='k', label='Species')

ax.legend(loc='upper left')

ax.set_xlabel('Mean relative error, synonymous SNVs', fontsize=11)
ax.set_ylabel('Mean relative error, nonsynonymous SNVs', fontsize=11)


fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%serror_4d_vs_1d.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
