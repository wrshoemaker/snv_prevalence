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
import scipy.stats as stats

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import calculate_predicted_prevalence_mapgd
import calculate_predicted_prevalence


import plot_utils



species_color_map, ordered_species_list = plot_utils.get_species_color_map()


prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()

clade_type = 'all'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'
max_n_occurances = 7


fig, ax_f = plt.subplots(figsize=(4,4))

for species_name in prevalence_dict_mapgd.keys():

    f_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_no_zeros_mapgd_folded']
    f_mapgd = numpy.asarray(f_mapgd)

    if len(f_mapgd) == 0:
        continue

    f_mapgd = 1-f_mapgd
    print(f_mapgd)

    print(species_name, sum(f_mapgd==1))

    f_mapgd_log10 = numpy.log10(f_mapgd)
    f_mapgd_log10_rescaled = (f_mapgd_log10 - numpy.mean(f_mapgd_log10)) / numpy.std(f_mapgd_log10)
    hist_f, bin_edges_f = numpy.histogram(f_mapgd_log10_rescaled, density=True, bins=20)
    bins_mean_f = [0.5 * (bin_edges_f[i] + bin_edges_f[i+1]) for i in range(0, len(bin_edges_f)-1 )]
    ax_f.scatter(bins_mean_f, hist_f, alpha=0.9, s=10, c=species_color_map[species_name])
    # label=figure_utils.get_pretty_species_name(species_name)


ax_f.set_xlabel('Rescaled ' + r'$\mathrm{log}_{10}$' + ' SNV frequencies', fontsize=11)
ax_f.set_ylabel('Probability density', fontsize=12)
#ax_f.set_ylim([-0.02, 1.23])



fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.25)
fig.savefig("%sdiversity_summary_folded_%s_%s.png" % (config.analysis_directory, clade_type, variant_type), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
