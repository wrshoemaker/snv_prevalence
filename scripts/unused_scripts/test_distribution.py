from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path
import random

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
import scipy.stats as stats
import scipy.special as special

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import calculate_predicted_prevalence_mapgd

import plot_utils



x = numpy.random.standard_normal(size=30)


x_gamma = numpy.random.gamma(7.5, scale=1.0, size=30)


fig, ax = plt.subplots(figsize=(4,4))

ax.hist(x_gamma, density=True, bins=20)


fig.tight_layout()
fig.subplots_adjust(wspace=0.34, hspace=0.28)
fig.savefig("%stest_dist.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()


n_draws = 1000
list_of_arrays = []

for n in range(n_draws):
    x_gamma_n = numpy.random.gamma(7.5, scale=1.0, size=30)
    list_of_arrays.append(x_gamma_n)

x_all = numpy.array(list_of_arrays).flatten()


fig, ax = plt.subplots(figsize=(4,4))

ax.hist(x_all, density=True, bins=20)

fig.tight_layout()
fig.subplots_adjust(wspace=0.34, hspace=0.28)
fig.savefig("%stest_dist_merged.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
