from __future__ import division
import sample_utils
import config
import parse_midas_data
import os.path
import os
import pylab
import sys
import numpy
import gzip
import pickle
import bz2

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
import prevalence_utils

from math import log10,ceil,fabs,isnan
from numpy.random import randint, choice, multinomial

import midas_db_utils

import scipy.stats as stats
from scipy.stats import gamma, gaussian_kde

import parse_HMP_data
import figure_utils
import calculate_predicted_occupancy

import matplotlib.pyplot as plt
import matplotlib.cm as cm


allowed_variant_types = set(['1D','4D'])

subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species = 'Eubacterium_rectale_56927'
#bad_species = 'Bacteroides_vulgatus_57955'


min_sample_size = config.between_host_min_sample_size
low_divergence_threshold = config.between_low_divergence_threshold

prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()


prevalence_min = 0.1
prevalence_max = 0.42
prevalence_range = numpy.logspace(-2.5, 0, num=1000)



fig = plt.figure(figsize = (8, 8)) #
fig.subplots_adjust(bottom= 0.15)

ax_prevalence_survival = plt.subplot2grid((2, 2), (0,0), colspan=2)

ax_good_species = plt.subplot2grid((2, 2), (1,0), colspan=1)

ax_good_species_fmax = plt.subplot2grid((2, 2), (1,1), colspan=1)


fontsize = 9

ax_prevalence_survival.text(0.01, 1.06, figure_utils.sub_plot_labels[0], fontsize=fontsize, fontweight='bold', ha='center', va='center', transform=ax_prevalence_survival.transAxes)

ax_good_species.text(0.01, 1.06, figure_utils.sub_plot_labels[1], fontsize=fontsize, fontweight='bold', ha='center', va='center', transform=ax_good_species.transAxes)
ax_good_species_fmax.text(0.01, 1.06, figure_utils.sub_plot_labels[2], fontsize=fontsize, fontweight='bold', ha='center', va='center', transform=ax_good_species_fmax.transAxes)


for species_name in prevalence_dict.keys():


    observed_prevalence = prevalence_dict[species_name]['all']['4D']['observed_prevalence_to_plot']
    observed_prevalence = numpy.asarray(observed_prevalence)

    survival_array = [sum(observed_prevalence>=i)/len(observed_prevalence) for i in prevalence_range]

    survival_array = numpy.asarray(survival_array)

    if species_name == good_species:
        color = prevalence_utils.good_bad_color_dict[species_name]
        zorder=2
        alpha=1
        lw=3.5
        ax_prevalence_survival.plot(prevalence_range, survival_array, ls='-', lw=lw, c=color, alpha=alpha, zorder=zorder, label = figure_utils.get_pretty_species_name(species_name))

    else:
        color = 'k'
        zorder=1
        alpha=0.25
        lw=2
        ax_prevalence_survival.plot(prevalence_range, survival_array, ls='-', lw=lw, c=color, alpha=alpha, zorder=zorder)





ax_prevalence_survival.legend(loc="lower left", fontsize=10)
ax_prevalence_survival.set_xlim([0.008, 1.05])
ax_prevalence_survival.set_xscale('log', basex=10)
ax_prevalence_survival.set_yscale('log', basey=10)
ax_prevalence_survival.set_xlabel('Observed SNV prevalence, ' + r'$\hat{\varrho}$', fontsize=14)
ax_prevalence_survival.set_ylabel('Fraction of SNVs ' + r'$\geq \hat{\varrho}$', fontsize=14)




for species_list in [[ax_good_species, ax_good_species_fmax, good_species]]:


    ax_species = species_list[0]
    ax_species_fmax = species_list[1]
    species_name = species_list[2]


    predicted_prevalence = prevalence_dict[species_name]['all']['4D']['predicted_prevalence_to_plot']
    predicted_prevalence = numpy.asarray(predicted_prevalence)

    observed_prevalence = prevalence_dict[species_name]['all']['4D']['observed_prevalence_to_plot']
    observed_prevalence = numpy.asarray(observed_prevalence)

    f_max = prevalence_dict[species_name]['all']['4D']['f_max']['to_plot']
    f_max = numpy.asarray(f_max)

    f_max_line = prevalence_dict[species_name]['all']['4D']['f_max_line']

    predicted_prevalence_line = prevalence_dict[species_name]['all']['4D']['f_max_vs_predicted_prevalence_line']

    f_max_line = numpy.asarray(f_max_line)
    predicted_prevalence_line = numpy.asarray(predicted_prevalence_line)


    # Calculate the point density
    xy = numpy.vstack([observed_prevalence,predicted_prevalence])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = observed_prevalence[idx], predicted_prevalence[idx], z[idx]


    ax_ = ax_species.scatter(x, y, c=z, cmap=prevalence_utils.good_bad_color_map_dict[species_name], s=70, alpha=0.9, edgecolor='')

    ax_species.scatter(10**-8, 10**-8, c=prevalence_utils.good_bad_color_dict[species_name], label='One SNV')


    ax_species.plot([min(observed_prevalence)*0.8, max(observed_prevalence)*1.1],[min(observed_prevalence)*0.8, max(observed_prevalence)*1.1], ls='--', lw=2, c='k', label='1:1')

    ax_species.set_xscale('log', basex=10)
    ax_species.set_yscale('log', basey=10)

    ax_species.set_xlabel('Observed SNV prevalence', fontsize=14)
    ax_species.set_ylabel('Predicted SNV prevalence', fontsize=14)


    ax_species.set_xlim([min(observed_prevalence)*0.8, max(observed_prevalence)*1.1])
    ax_species.set_ylim([min(observed_prevalence)*0.8, max(observed_prevalence)*1.1])


    ax_species.set_title(figure_utils.get_pretty_species_name(species_name), fontsize=12, fontweight='bold', color='k' )

    ax_species.legend(loc="upper left", fontsize=10)




    all_ = numpy.concatenate([f_max,observed_prevalence])


    # Calculate the point density
    xy = numpy.vstack([f_max, observed_prevalence])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = f_max[idx], observed_prevalence[idx], z[idx]

    #ax.set_xlim([min(f_max)*0.8, max(f_max)*1.1])
    ax_species_fmax.set_xlim([min(f_max)*0.8, 0.2])
    ax_species_fmax.set_ylim([min(observed_prevalence)*0.8, max(observed_prevalence)*1.1])

    ax_species_fmax.scatter(x, y, c=z, cmap=prevalence_utils.good_bad_color_map_dict[species_name], s=25, alpha=0.9, edgecolor='', zorder=1)

    ax_species_fmax.plot(10**f_max_line, 10**predicted_prevalence_line, c='k', lw = 3, ls='--', label = 'Gamma prediction', zorder=2)
    ax_species_fmax.scatter(10**-5, 10**-5, c=prevalence_utils.good_bad_color_dict[species_name], label='Observation of one SNV')

    ax_species_fmax.set_title(figure_utils.get_pretty_species_name(species_name), fontsize=12, fontweight='bold', color='k' )

    ax_species_fmax.set_xscale('log', basex=10)
    ax_species_fmax.set_yscale('log', basey=10)

    ax_species_fmax.set_xlabel('Maximum frequency\nacross hosts, ' + r'$f_{max}$', fontsize=14)
    ax_species_fmax.set_ylabel('SNV prevalence', fontsize=14)


    ax_species_fmax.legend(loc="upper left", fontsize=10)



fig.tight_layout()
fig.savefig("%sprevalence_summary.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
