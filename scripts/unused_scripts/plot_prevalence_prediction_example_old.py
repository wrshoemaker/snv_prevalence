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
import scipy.special as special

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_prevalence_mapgd

import plot_utils


prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all(test=True)
species_color_map, ordered_species_list = plot_utils.get_species_color_map()

clade_type = 'all'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
variant_type = '4D'
species_name = 'Bacteroides_vulgatus_57955'

fig = plt.figure(figsize = (12, 8)) #

ax_prevalence_eco = plt.subplot2grid((2, 3), (0, 0), colspan=1)
ax_prevalence_evo = plt.subplot2grid((2, 3), (1, 0), colspan=1)
ax_abundance_prevalence_eco = plt.subplot2grid((2, 3), (0, 1), colspan=1)
ax_abundance_prevalence_evo = plt.subplot2grid((2, 3), (1, 1), colspan=1)
ax_error = plt.subplot2grid((2, 3), (0, 2), colspan=1)
ax_error_vs = plt.subplot2grid((2, 3), (1, 2), colspan=1)

ax_all = [ax_prevalence_eco, ax_prevalence_evo, ax_abundance_prevalence_eco, ax_abundance_prevalence_evo, ax_error, ax_error_vs]
for ax_idx, ax in enumerate(ax_all):
    ax.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[ax_idx], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax.transAxes)



observed_prevalence_slm = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
observed_prevalence_slm = numpy.asarray(observed_prevalence_slm)

predicted_prevalence_slm = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm']
predicted_prevalence_slm = numpy.asarray(predicted_prevalence_slm)

observed_prevalence = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
observed_prevalence = numpy.asarray(observed_prevalence)

predicted_prevalence = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd']
predicted_prevalence = numpy.asarray(predicted_prevalence)

f_mean = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_mean_mapgd']
f_mean = numpy.asarray(f_mean)

f_max = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_max_mapgd']
f_max = numpy.asarray(f_max)



idx_ = (observed_prevalence>0) & (predicted_prevalence>0) & (observed_prevalence_slm>0) & (predicted_prevalence_slm>0) & (f_mean>0) & (f_max>0)

observed_prevalence_slm = observed_prevalence_slm[idx_]
predicted_prevalence_slm = predicted_prevalence_slm[idx_]
observed_prevalence = observed_prevalence[idx_]
predicted_prevalence = predicted_prevalence[idx_]
f_mean = f_mean[idx_]
f_max = f_max[idx_]


# eco
all_ = numpy.concatenate([predicted_prevalence_slm,observed_prevalence_slm])
xy = numpy.vstack([observed_prevalence_slm, predicted_prevalence_slm])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last
idx_ = z.argsort()
x, y, z = observed_prevalence_slm[idx_], predicted_prevalence_slm[idx_], z[idx_]
ax_prevalence_eco.scatter(x, y, c=z, cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolor='', zorder=1)

#sorted_plot_data = prevalence_utils.plot_color_by_pt_dens(observed_prevalence_slm, predicted_prevalence_slm, radius=prevalence_utils.color_radius, loglog=1)
#x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]
#ax_prevalence_eco.scatter(x, y, c=numpy.sqrt(z), cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolors='none', zorder=1)
#all_ = numpy.concatenate([x, y])

max_ = max(all_)*1.2
min_ = min(all_)*0.8

ax_prevalence_eco.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2, label='1:1')
ax_prevalence_eco.set_xlim([min_, max_])
ax_prevalence_eco.set_ylim([min_, max_])
ax_prevalence_eco.set_xscale('log', basex=10)
ax_prevalence_eco.set_yscale('log', basey=10)
ax_prevalence_eco.set_title("Stochastic Logistic Model", fontsize=12, fontweight='bold', color='k' )
ax_prevalence_eco.set_xlabel('Observed allele prevalence', fontsize=12)
ax_prevalence_eco.set_ylabel('Predicted allele prevalence', fontsize=12)
ax_prevalence_eco.legend(loc="upper left", fontsize=12)
ax_prevalence_eco.xaxis.set_tick_params(labelsize=8)
ax_prevalence_eco.yaxis.set_tick_params(labelsize=8)



# evo
all_ = numpy.concatenate([predicted_prevalence,observed_prevalence])
xy = numpy.vstack([observed_prevalence, predicted_prevalence])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last
idx_ = z.argsort()
x, y, z = observed_prevalence[idx_], predicted_prevalence[idx_], z[idx_]
ax_prevalence_evo.scatter(x, y, c=z, cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolor='', zorder=1)

#sorted_plot_data = prevalence_utils.plot_color_by_pt_dens(observed_prevalence, predicted_prevalence, radius=prevalence_utils.color_radius, loglog=1)
#x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]
#ax_prevalence_evo.scatter(x, y, c=numpy.sqrt(z), cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolors='none', zorder=1)
#all_ = numpy.concatenate([x, y])



max_ = max(all_)*1.2
min_ = min(all_)*0.8

ax_prevalence_evo.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2, label='1:1')
ax_prevalence_evo.set_xlim([min_, max_])
ax_prevalence_evo.set_ylim([min_, max_])
ax_prevalence_evo.set_xscale('log', basex=10)
ax_prevalence_evo.set_yscale('log', basey=10)
ax_prevalence_evo.set_title("Evolutionary model", fontsize=12, fontweight='bold', color='k' )
ax_prevalence_evo.set_xlabel('Observed allele prevalence', fontsize=12)
ax_prevalence_evo.set_ylabel('Predicted allele prevalence', fontsize=12)
ax_prevalence_evo.legend(loc="upper left", fontsize=12)
ax_prevalence_evo.xaxis.set_tick_params(labelsize=8)
ax_prevalence_evo.yaxis.set_tick_params(labelsize=8)



# eco, abundnace occupancy
f_mean_log10 = numpy.log10(f_mean)
# relationship between f_mean and prevalence
hist, bin_edges = numpy.histogram(f_mean_log10, density=True, bins=40)
bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]
prevalence_predicted_mean = []
for i in range(0, len(bin_edges)-1):
    idx_i = (f_mean_log10 > bin_edges[i]) & (f_mean_log10 <= bin_edges[i+1])
    prevalence_predicted_mean.append(numpy.mean(numpy.log10(predicted_prevalence_slm[idx_i])))
bins_mean = numpy.asarray(bins_mean)
prevalence_predicted_mean = numpy.asarray(prevalence_predicted_mean)



all_ = numpy.concatenate([f_mean, observed_prevalence_slm])
xy = numpy.vstack([f_mean, observed_prevalence_slm])
z = gaussian_kde(xy)(xy)
#max_ = max(all_)*1.2
#min_ = min(all_)*0.8
#min_ = 0.003864734299516908
# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = f_mean[idx], observed_prevalence_slm[idx], z[idx]
ax_abundance_prevalence_eco.scatter(x, y, c=z, cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=1, edgecolor='', zorder=1, label = 'Observed')

#sorted_plot_data = prevalence_utils.plot_color_by_pt_dens(f_mean, observed_prevalence_slm, radius=prevalence_utils.color_radius, loglog=1)
#x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]
#ax_abundance_prevalence_eco.scatter(x, y, c=numpy.sqrt(z), cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolors='none', zorder=1, label = 'Observed')
#all_ = numpy.concatenate([x, y])

max_ = max(all_)*1.2
#min_ = min(all_)*0.8
min_ = 0.003864734299516908

ax_abundance_prevalence_eco.plot(10**bins_mean, 10**prevalence_predicted_mean, c='k', lw = 3, ls='--', label = 'Predicted', zorder=2)
ax_abundance_prevalence_eco.set_title("Stochastic Logistic Model", fontsize=12, fontweight='bold', color='k' )
ax_abundance_prevalence_eco.set_xlim([min_, max_])
ax_abundance_prevalence_eco.set_ylim([min_, max_])
ax_abundance_prevalence_eco.set_xscale('log', basex=10)
ax_abundance_prevalence_eco.set_yscale('log', basey=10)
ax_abundance_prevalence_eco.set_xlabel('Mean within-host allele\nfrequency across hosts, ' + r'$\bar{f}$', fontsize=12)
ax_abundance_prevalence_eco.set_ylabel('Allele prevalence', fontsize=12)
ax_abundance_prevalence_eco.legend(loc="upper left", fontsize=12)
ax_abundance_prevalence_eco.xaxis.set_tick_params(labelsize=8)
ax_abundance_prevalence_eco.yaxis.set_tick_params(labelsize=8)




# evo, abundnace occupancy
f_max_log10 = numpy.log10(f_max)
# relationship between f_max and prevalence
hist, bin_edges = numpy.histogram(f_max_log10, density=True, bins=40)
bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]
prevalence_predicted_mean = []
for i in range(0, len(bin_edges)-1):
    idx_i = (f_max_log10 > bin_edges[i]) & (f_max_log10 <= bin_edges[i+1])
    prevalence_predicted_mean.append(numpy.mean(numpy.log10(predicted_prevalence[idx_i])))
bins_mean = numpy.asarray(bins_mean)
prevalence_predicted_mean = numpy.asarray(prevalence_predicted_mean)

all_ = numpy.concatenate([f_max, observed_prevalence])
xy = numpy.vstack([f_max, observed_prevalence])
z = gaussian_kde(xy)(xy)

#sorted_plot_data = prevalence_utils.plot_color_by_pt_dens(f_max, observed_prevalence, radius=prevalence_utils.color_radius, loglog=1)
#x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]
#ax_abundance_prevalence_evo.scatter(x, y, c=numpy.sqrt(z), cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolors='none', zorder=1, label='Observed')
#all_ = numpy.concatenate([x, y])


max_ = max(all_)*1.2
min_ = min(all_)*0.8
# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = f_max[idx], observed_prevalence[idx], z[idx]
ax_abundance_prevalence_evo.scatter(x, y, c=z, cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=1, edgecolor='', zorder=1, label='Observed')
ax_abundance_prevalence_evo.plot(10**bins_mean, 10**prevalence_predicted_mean, c='k', lw = 3, ls='--', label = 'Predicted', zorder=2)
ax_abundance_prevalence_evo.set_title("Evolutionary Model", fontsize=12, fontweight='bold', color='k' )
ax_abundance_prevalence_evo.set_xlim([min_, max_])
ax_abundance_prevalence_evo.set_ylim([min_, max_])
ax_abundance_prevalence_evo.set_xscale('log', basex=10)
ax_abundance_prevalence_evo.set_yscale('log', basey=10)
ax_abundance_prevalence_evo.set_xlabel('Max. within-host allele\nfrequency across hosts, ' + r'$f_{max}$', fontsize=12)
ax_abundance_prevalence_evo.set_ylabel('Allele prevalence', fontsize=12)
ax_abundance_prevalence_evo.legend(loc="upper left", fontsize=12)
ax_abundance_prevalence_evo.xaxis.set_tick_params(labelsize=8)
ax_abundance_prevalence_evo.yaxis.set_tick_params(labelsize=8)


# error survival distribution
error_eco = numpy.absolute(predicted_prevalence_slm - observed_prevalence_slm) / observed_prevalence_slm
error_evo = numpy.absolute(predicted_prevalence - observed_prevalence) / observed_prevalence

error_eco_log10 = numpy.log10(error_eco)
error_evo_log10 = numpy.log10(error_evo)

idx_ = (numpy.isfinite(error_eco_log10)) & (numpy.isfinite(error_evo_log10))

error_eco_log10_no_nan = error_eco_log10[idx_]
error_evo_log10_no_nan = error_evo_log10[idx_]


error_eco_no_nan = 10**error_eco_log10_no_nan
error_evo_no_nan = 10**error_evo_log10_no_nan

min_x = min([min(error_eco_no_nan), min(error_evo_no_nan)])
max_x = max([max(error_eco_no_nan), max(error_evo_no_nan)])

x_range = numpy.logspace(numpy.log10(min_x), numpy.log10(max_x), num=100, endpoint=True, base=10.0)

survival_error_eco = [sum(error_eco_no_nan >= i)/len(error_eco_no_nan) for i in x_range]
survival_error_evo = [sum(error_evo_no_nan >= i)/len(error_evo_no_nan) for i in x_range]

survival_error_eco = numpy.asarray(survival_error_eco)
survival_error_evo = numpy.asarray(survival_error_evo)


ax_error.axvline(10**numpy.mean(error_eco_log10_no_nan), lw=1.5, ls='--',color='k', zorder=1, label='SLM, ' + r'$\overline{\mathrm{log}_{10} \epsilon}$')
ax_error.axvline(10**numpy.mean(error_evo_log10_no_nan), lw=1.5, ls=':',color='k', zorder=1, label='Evo model, ' + r'$\overline{\mathrm{log}_{10} \epsilon}$')

ax_error.plot(x_range, survival_error_eco, ls='--', lw=2, c=species_color_map[species_name], zorder=2, label='SLM')
ax_error.plot(x_range, survival_error_evo, ls=':', lw=2, c=species_color_map[species_name], zorder=2, label='Evo model')
ax_error.set_xscale('log', basex=10)
ax_error.legend(loc="lower left", fontsize=8)
ax_error.set_xlabel('Relative error of prevalence prediction, ' + r'$\epsilon$', fontsize=12)
ax_error.set_ylabel('Fraction of sites ' + r'$\geq \epsilon$', fontsize=12)
ax_error.xaxis.set_tick_params(labelsize=8)
ax_error.yaxis.set_tick_params(labelsize=8)


# error comparison
all_ = numpy.concatenate([error_evo_no_nan,error_eco_no_nan])
xy = numpy.vstack([error_eco_no_nan, error_evo_no_nan])
z = gaussian_kde(xy)(xy)
# Sort the points by density, so that the densest points are plotted last
idx_ = z.argsort()
x, y, z = error_eco_no_nan[idx_], error_evo_no_nan[idx_], z[idx_]
ax_error_vs.scatter(x, y, c=z, cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolor='', zorder=1)

#sorted_plot_data = prevalence_utils.plot_color_by_pt_dens(error_eco_no_nan, error_evo_no_nan, radius=prevalence_utils.color_radius, loglog=1)
#x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]
#ax_error_vs.scatter(x, y, c=numpy.sqrt(z), cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolors='none', zorder=1)
#all_ = numpy.concatenate([x, y])

max_ = max(all_)*1.2
min_ = min(all_)*0.8

ax_error_vs.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2, label='1:1')
ax_error_vs.set_xlim([min_, max_])
ax_error_vs.set_ylim([min_, max_])
ax_error_vs.set_xscale('log', basex=10)
ax_error_vs.set_yscale('log', basey=10)
ax_error_vs.set_xlabel('Relative error, SLM', fontsize=12)
ax_error_vs.set_ylabel('Relative error, evo model', fontsize=12)
ax_error_vs.legend(loc="lower right", fontsize=12)
ax_error_vs.xaxis.set_tick_params(labelsize=8)
ax_error_vs.yaxis.set_tick_params(labelsize=8)





fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.28)
fig.savefig("%sprevalence_example.pdf" % config.analysis_directory, format='pdf', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
