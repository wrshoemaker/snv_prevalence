import matplotlib
matplotlib.use('Agg')
import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
from math import exp

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_temporal_changes
import calculate_substitution_rates
import calculate_linkage_disequilibria
import calculate_snv_distances
import stats_utils
import sfs_utils
import figure_utils

import scipy.stats as stats


from scipy.optimize import least_squares, newton, brentq

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, multinomial
import matplotlib.colors as mcolors

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import calculate_predicted_occupancy

from math import log


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()

species_list = list(prevalence_dict.keys())
species_list.sort()


fig = plt.figure(figsize = (12, 4)) #
fig.subplots_adjust(bottom= 0.15)


ax_9 = plt.subplot2grid((1, 3), (0,0))
ax_90 = plt.subplot2grid((1, 3), (0,1))
ax_900 = plt.subplot2grid((1, 3), (0,2))

cmap_offset_co = int(0.4*1000)
rgb_blue_co = cm.Blues(numpy.linspace(0,1,1000+cmap_offset_co ))
rgb_blue_co = mpl.colors.ListedColormap(rgb_blue_co[cmap_offset_co:,:-1])


ld_plot_dict = {}
ld_plot_dict['mre'] = []
ld_plot_dict['9']  = []
ld_plot_dict['90']  = []
ld_plot_dict['900']  = []

for species_idx, species_name in enumerate(species_list):
    # Load precomputed LD
    ld_map = calculate_linkage_disequilibria.load_ld_map(species_name)

    if len(ld_map)==0:
        continue

    all_distances, all_rsquared_numerators, all_rsquared_denominators, all_ns, all_intergene_distances, all_intergene_rsquared_numerators, all_intergene_rsquared_denominators, all_intergene_ns, all_control_rsquared_numerator, all_control_rsquared_denominator, all_control_n, all_pi = ld_map[('all','4D')]
    all_control_rsquared = all_control_rsquared_numerator/all_control_rsquared_denominator

    distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerator, control_rsquared_denominator, control_n, pi = ld_map[('largest_clade','4D')]
    control_rsquared = control_rsquared_numerator/control_rsquared_denominator

    # smooth this stuff:
    smoothed_distances = distances
    window_width = 10**(0.1)

    dmins = smoothed_distances/(window_width**0.5)
    dmaxs = smoothed_distances*(window_width**0.5)

    smoothed_rsquared_numerators = []
    smoothed_rsquared_denominators = []
    smoothed_counts = []

    all_smoothed_rsquared_numerators = []
    all_smoothed_rsquared_denominators = []
    all_smoothed_counts = []

    for dmin,dmax in zip(dmins,dmaxs):
        binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
        binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
        binned_counts = ns[(distances>=dmin)*(distances<=dmax)]
        smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
        smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
        smoothed_counts.append( binned_counts.sum() )

        binned_numerators = all_rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
        binned_denominators = all_rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
        binned_counts = all_ns[(distances>=dmin)*(distances<=dmax)]
        all_smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
        all_smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
        all_smoothed_counts.append( binned_counts.sum() )


    smoothed_rsquared_numerators = numpy.array( smoothed_rsquared_numerators )
    smoothed_rsquared_denominators = numpy.array( smoothed_rsquared_denominators )
    smoothed_counts = numpy.array( smoothed_counts )

    all_smoothed_rsquared_numerators = numpy.array( all_smoothed_rsquared_numerators )
    all_smoothed_rsquared_denominators = numpy.array( all_smoothed_rsquared_denominators )
    all_smoothed_counts = numpy.array( all_smoothed_counts )

    early_distances = distances[distances<101]
    early_rsquareds = rsquared_numerators[distances<101]*1.0/rsquared_denominators[distances<101]
    early_ns = ns[distances<101]

    early_distances = early_distances[early_ns>0.5]
    early_rsquareds = early_rsquareds[early_ns>0.5]
    early_ns = early_ns[early_ns>0.5]

    distances = smoothed_distances
    rsquareds = smoothed_rsquared_numerators/(smoothed_rsquared_denominators)
    ns = smoothed_counts

    distances = distances[ns>0]
    rsquareds = rsquareds[ns>0]
    ns = ns[ns>0]

    all_distances = smoothed_distances
    #all_distances = dmins
    all_rsquareds = all_smoothed_rsquared_numerators/(all_smoothed_rsquared_denominators)
    all_ns = all_smoothed_counts

    all_distances = all_distances[all_ns>0]
    all_rsquareds = all_rsquareds[all_ns>0]
    all_ns = all_ns[all_ns>0]


    idx_9 = numpy.fabs(distances-9).argmin()
    idx_90 = numpy.fabs(distances-100).argmin()
    idx_900 = numpy.fabs(distances-2000).argmin()

    rsquared_9 = rsquareds[idx_9]
    rsquared_90 = rsquareds[idx_90]
    rsquared_900 = rsquareds[idx_900]


    mre = prevalence_dict[species_name]['all']['4D']['MRE']

    ld_plot_dict['mre'].append(mre)
    ld_plot_dict['9'].append(rsquared_9)
    ld_plot_dict['90'].append(rsquared_90)
    ld_plot_dict['900'].append(rsquared_900)

    ax_9.scatter(rsquared_9, mre, c=rgb_blue_co(idx_9), edgecolors='k', alpha=1)
    ax_90.scatter(rsquared_90, mre, c=rgb_blue_co(idx_90), edgecolors='k', alpha=1)
    ax_900.scatter(rsquared_900, mre, c=rgb_blue_co(idx_900), edgecolors='k', alpha=1)




mre = ld_plot_dict['mre']
ld_9 = ld_plot_dict['9']
ld_90 = ld_plot_dict['90']
ld_900 = ld_plot_dict['900']


mre = numpy.asarray(mre)
ld_9 = numpy.asarray(ld_9)
ld_90 = numpy.asarray(ld_90)
ld_900 = numpy.asarray(ld_900)



slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(ld_9), mre)
ax_9.text(0.85,0.9, r'$P \nless 0.05$', fontsize=12, color='k', ha='center', va='center', transform=ax_9.transAxes  )


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(ld_90), mre)
ax_90.text(0.85,0.9, r'$P \nless 0.05$', fontsize=12, color='k', ha='center', va='center', transform=ax_90.transAxes  )

slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(ld_900), mre)
ax_900.text(0.85,0.9, r'$P \nless 0.05$', fontsize=12, color='k', ha='center', va='center', transform=ax_900.transAxes  )


ax_9.set_xlabel('Linkage disequilibrium, $\sigma^2_d$')
ax_90.set_xlabel('Linkage disequilibrium, $\sigma^2_d$')
ax_900.set_xlabel('Linkage disequilibrium, $\sigma^2_d$')

ax_9.set_ylabel('Mean relative error')
ax_90.set_ylabel('Mean relative error')
ax_900.set_ylabel('Mean relative error')


ax_9.set_xscale('log', basex=10)
ax_90.set_xscale('log', basex=10)
ax_900.set_xscale('log', basex=10)


ax_9.set_title(r'$\Delta \ell = 10$', fontsize=12, fontweight='bold', color='k' )
ax_90.set_title(r'$\Delta \ell = 100$', fontsize=12, fontweight='bold', color='k' )
ax_900.set_title(r'$\Delta \ell = 1,000$', fontsize=12, fontweight='bold', color='k' )



fig.subplots_adjust(hspace=0.4, wspace=0.4)
fig.savefig("%smre_vs_ld.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
