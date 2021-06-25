import sys, os
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria

import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#import calculate_snv_distances
import figure_utils
from math import log10,ceil
#import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.random import randint, multinomial
import parse_HMP_data

import figure_utils
import ld_utils

import matplotlib as mpl


bin_width_exponent = 1.1
variant_types = ['4D','1D']

import statsmodels.api as sm


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize



sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species_list = parse_midas_data.parse_good_species_list()

ld_dict_all_species = calculate_linkage_disequilibria.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)

ld_species = list(ld_dict_all_species.keys())

n_cols = 6
species_nested_list  = list(ld_utils.chunks(ld_species, n_cols))
n_rows = len(species_nested_list)



fig = plt.figure(figsize = (10, 8))
gs = gridspec.GridSpec(nrows=n_rows, ncols=n_cols)

color_dict = {'4D': '#87CEEB', '1D':'#FF6347'}

x_axis_label = 'Change in LD ratio\nover $\Delta \ell$bp, $\Delta \sigma^2_{d,N} / \sigma^2_{d,S}$'


fig, ax = plt.subplots(figsize=(4,4))

dmin=10
dmax=1000

all_x = []
all_y = []


ax.plot([0.001,100],[0.001,100], lw=2,ls='--',c='k',zorder=1)

ax.plot([0.001,1],[1,1], lw=2,ls=':',c='k',zorder=1)
ax.plot([1,1],[0.001,1], lw=2,ls=':',c='k',zorder=1)
ax.set_xlabel("Mean LD ratio across bins", fontsize = 12)
ax.set_ylabel("Genome-wide LD ratio from clones", fontsize = 12)

for species in ld_species:

    ld_dict = ld_dict_all_species[species]

    ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, control_rsquared_ratio = ld_utils.calculate_ld_4D_1D_ratio(ld_dict)

    ld_1D_4D_ratio_filter_rsquareds_dmin_dmax = ld_1D_4D_ratio_filter_rsquareds[(ld_filter_distances>dmin) & (ld_filter_distances<dmax)]

    mean_ld_1D_4D_ratio_filter_rsquareds = numpy.mean(ld_1D_4D_ratio_filter_rsquareds_dmin_dmax)

    all_x.append(mean_ld_1D_4D_ratio_filter_rsquareds)
    all_y.append(control_rsquared_ratio)
    ax.scatter(mean_ld_1D_4D_ratio_filter_rsquareds, control_rsquared_ratio, zorder=2)



all_x = numpy.asarray(all_x)
all_y = numpy.asarray(all_y)

max_x_idx = numpy.where(all_x == numpy.amax(all_x))
min_x_idx = numpy.where(all_x == numpy.amin(all_x))
max_species = ld_species[max_x_idx[0][0]]
min_species = ld_species[min_x_idx[0][0]]

print('max ', max_species)
print('min ', min_species)

#model = sm.OLS(all_x, all_y)
#results = model.fit()
#print(results.summary())


#plt.hlines(y=1, xmin=0.001, xmax=0, color='k', linestyle=':', lw = 2, zorder=1)
#plt.vlines(x=1, ymin=0.001, ymax=0, color='k', linestyle=':', lw = 2, zorder=1)


max_lim = max(max(all_x), max(all_y))
min_lim = min(min(all_x), min(all_y))

ax.set_xlim([min_lim*0.8,max_lim*1.1])
ax.set_ylim([min_lim*0.8,max_lim*1.1])


#ax.set_xscale('log', basex=10)
#ax.set_yscale('log', basey=10)

# plot of the LD distance decay relationships for all species

fig.subplots_adjust(hspace=0.15, wspace=0.15)
fig.savefig("%s%s.png" % (config.analysis_directory, 'ld_ratio_vs_control'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
