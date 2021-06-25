import sys, os, re
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria
from math import log10,ceil

#import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

#import calculate_snv_distances
import figure_utils
from math import log10,ceil
from numpy.random import randint, multinomial
import parse_HMP_data

import figure_utils
import ld_utils
import matplotlib as mpl
from matplotlib import cm

import scipy.stats as stats


import calculate_predicted_occupancy

prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()

#good_species_list = prevalence_dict.keys()

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

ld_dict_all_species = calculate_linkage_disequilibria.calculate_ld_dict(prevalence_dict.keys(), subject_sample_map)

#print(ld_dict_all_species)

good_species_list = list(set(prevalence_dict.keys()) & set(ld_dict_all_species.keys()) )


mae_ld_all = []
mae_prevalence_all = []

fig, ax = plt.subplots(figsize=(4,4))

for species_name in good_species_list:

    ld_dict = ld_dict_all_species[species_name]

    distances = ld_dict['4D']['distances']
    rsquareds = ld_dict['4D']['rsquareds']

    upper_rsquareds = ld_dict['4D']['upper_rsquareds']
    lower_rsquareds = ld_dict['4D']['lower_rsquareds']


    good_distances = (upper_rsquareds>=-0.5)*(lower_rsquareds>=-0.5)
    #theory_ls = numpy.logspace(0,log10(distances[-1]),100)
    theory_ls = distances[good_distances]
    theory_NRs = theory_ls/200.0
    theory_rsquareds = (10+2*theory_NRs)/(22+26*theory_NRs+4*theory_NRs*theory_NRs)

    theory_rsquareds = (theory_rsquareds/theory_rsquareds[0]*3e-0)
    #print(theory_rsquareds)

    rsquareds = rsquareds[good_distances]

    distances = distances[good_distances]

    mae_ld = numpy.mean(numpy.absolute(rsquareds - theory_rsquareds) / distances)

    #mae_ld = numpy.mean(rsquareds/distances)

    mae_prevalence = prevalence_dict[species_name]['predicted_observed_prevalence']['MAE']

    ax.scatter(mae_prevalence, mae_ld)

    mae_ld_all.append(mae_ld)
    mae_prevalence_all.append(mae_prevalence)



mae_ld_all = numpy.asarray(mae_ld_all)
mae_prevalence_all = numpy.asarray(mae_prevalence_all)


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mae_prevalence_all), numpy.log10(mae_ld_all))

print(p_value)


    #example_axis.loglog(theory_ls, theory_rsquareds/theory_rsquareds[0]*3e-01,'k-',linewidth=0.3,zorder=0,label='Neutral')

ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel('Mean absolute error of\npredicted/observed prevalence')
ax.set_ylabel('Mean absolute error of\npredicted/observed per-bp LD')

ax.set_ylim([0.003, 0.006 ])
ax.set_xlim([0.01, 0.13 ])

#ax.xaxis.set_tick_params(labelsize=6)
#ax.yaxis.set_tick_params(labelsize=6)

#ax.tick_params(axis='both', which='major', labelsize=10)
ax.tick_params(axis='both', which='minor', labelsize=7)
ax.tick_params(axis='both', which='major', labelsize=7)


fig.tight_layout()
fig.savefig("%sprevalence_error_vs_ld_error.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
