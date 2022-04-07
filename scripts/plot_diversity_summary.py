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

import calculate_predicted_prevalence_mapgd

import plot_utils

#species_color_map, ordered_species_list = plot_utils.get_species_color_map()
species_color_map = prevalence_utils.species_color_map_genus

prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()#test=True)

species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()

clade_type = 'all'
#clade_type = 'all'
pi_type = 'pi_include_boundary'
f_mean_upper_cutoff_taylor = 0.35
max_n_occurances = 10
max_n_occurances_range = list(range(1, max_n_occurances+1))




def get_f_no_zeros_dict(species_name, clade_type='all'):

    intermediate_filename_template = config.data_directory+"frequency_dict_mapgd_non_zero/%s_%s.dat"
    intermediate_filename = intermediate_filename_template % (species_name, clade_type)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b


def make_plot(variant_type):

    fig = plt.figure(figsize = (12, 8)) #

    ax_n_hosts = plt.subplot2grid((2, 3), (0, 0), colspan=1)
    ax_n_sites = plt.subplot2grid((2, 3), (1, 0), colspan=1)
    ax_f = plt.subplot2grid((2, 3), (0, 1), colspan=1)
    ax_f_mean = plt.subplot2grid((2, 3), (1, 1), colspan=1)
    ax_f_mean_vs_var = plt.subplot2grid((2, 3), (0, 2), colspan=1)
    ax_f_prevalence = plt.subplot2grid((2, 3), (1, 2), colspan=1)

    ax_all = [ax_n_hosts, ax_n_sites, ax_f, ax_f_mean, ax_f_mean_vs_var, ax_f_prevalence]

    for ax_idx, ax in enumerate(ax_all):
        ax.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[ax_idx], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax.transAxes)


    if clade_type == 'all':
        min_prevalence = 0.35
    else:
        min_prevalence = 0.1


    relative_richness_dict, n_samples_dict, proprtion_richness = prevalence_utils.get_relative_richness_dict(variant_type=variant_type)


    # plot number of hosts
    species_n_hosts_list = [(species_name, n_samples_dict[species_name]) for species_name in species_list]
    species_n_hosts_list.sort(key=lambda tup: tup[1])
    species_n_hosts_list = species_n_hosts_list[::-1]
    species_list_hosts_to_plot = [s[0] for s in species_n_hosts_list]
    hosts_to_plot = [s[1] for s in species_n_hosts_list]

    species_list_hosts_to_plot_pretty = [figure_utils.get_pretty_species_name(s) for s in species_list_hosts_to_plot]
    colors_hosts = [species_color_map[s] for s in species_list_hosts_to_plot]
    ax_n_hosts.barh(species_list_hosts_to_plot_pretty, hosts_to_plot, height=0.8, align='center', color=colors_hosts)
    ax_n_hosts.set_xlabel('Number of hosts', fontsize=12)
    ax_n_hosts.xaxis.set_tick_params(labelsize=8)
    ax_n_hosts.yaxis.set_tick_params(labelsize=7)
    ax_n_hosts.set_ylim([-0.6, len(species_list_hosts_to_plot_pretty)-0.3])



    sys.stderr.write("%s sites\n" % variant_type)

    sys.stderr.write("Min. # hosts = %d\n" % min(hosts_to_plot))
    sys.stderr.write("Max. # hosts = %d\n" % max(hosts_to_plot))
    sys.stderr.write("Median # hosts = %f\n" % numpy.median(hosts_to_plot))

    # plot number of snvs
    # sort species names
    species_n_sites_list = []
    for species_name in species_list:
        observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
        observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)
        #n_sites = sum(observed_prevalence_mapgd >= min_prevalence)
        n_sites = len(observed_prevalence_mapgd)
        species_n_sites_list.append((species_name, n_sites))

    species_n_sites_list.sort(key=lambda tup: tup[1])
    species_n_sites_list = species_n_sites_list[::-1]
    species_list_sites_to_plot = [s[0] for s in species_n_sites_list]
    n_sites_to_plot = [s[1] for s in species_n_sites_list]

    species_list_sites_to_plot_pretty = [figure_utils.get_pretty_species_name(s) for s in species_list_sites_to_plot]
    colors_sites = [species_color_map[s] for s in species_list_hosts_to_plot]
    ax_n_sites.barh(species_list_sites_to_plot_pretty, n_sites_to_plot, height=0.8, align='center', color=colors_sites)
    ax_n_sites.set_xlabel('Number of sites', fontsize=12)
    ax_n_sites.xaxis.set_tick_params(labelsize=8)
    ax_n_sites.yaxis.set_tick_params(labelsize=7)
    ax_n_sites.set_ylim([-0.6, len(species_list_sites_to_plot_pretty)-0.3])
    ax_n_sites.set_xscale('log', basex =10)

    sys.stderr.write("Min. # sites = %d\n" % min(n_sites_to_plot))
    sys.stderr.write("Max. # sites = %d\n" % max(n_sites_to_plot))
    sys.stderr.write("Median # sites = %f\n" % 10**numpy.median(numpy.log10(n_sites_to_plot)))


    # plot everything else

    means_all = []
    variances_all = []
    prevalence_all = []
    n_non_zero_f_all = []
    means_all_for_prevalence = []
    species_list_to_plot_good = []

    min_all = []
    for species_name in species_list_sites_to_plot:

        f_no_zeros_mapgd_dict = get_f_no_zeros_dict(species_name)
        f_no_zeros_mapgd = f_no_zeros_mapgd_dict[variant_type]
        f_no_zeros_mapgd = numpy.asarray(f_no_zeros_mapgd)

        min_all.append(min(f_no_zeros_mapgd))

        n_non_zero_f = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
        n_non_zero_f = numpy.asarray(n_non_zero_f)

        f_mean_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_mean_mapgd']
        f_mean_mapgd = numpy.asarray(f_mean_mapgd)

        f_var_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_var_mapgd']
        f_var_mapgd = numpy.asarray(f_var_mapgd)

        observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
        observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)

        if len(f_no_zeros_mapgd) == 0:
            continue

        f_mean_mapgd_to_plot = f_mean_mapgd[observed_prevalence_mapgd >= min_prevalence]
        f_var_mapgd_to_plot = f_var_mapgd[observed_prevalence_mapgd >= min_prevalence]

        f_no_zeros_mapgd_log10 = numpy.log10(f_no_zeros_mapgd)
        f_no_zeros_mapgd_log10_rescaled = (f_no_zeros_mapgd_log10 - numpy.mean(f_no_zeros_mapgd_log10)) / numpy.std(f_no_zeros_mapgd_log10)
        hist_f, bin_edges_f = numpy.histogram(f_no_zeros_mapgd_log10_rescaled, density=True, bins=20)
        bins_mean_f = [0.5 * (bin_edges_f[i] + bin_edges_f[i+1]) for i in range(0, len(bin_edges_f)-1 )]
        bins_mean_f = numpy.asarray(bins_mean_f)
        hist_f_to_plot = hist_f[hist_f>0]
        bins_mean_f_to_plot = bins_mean_f[hist_f>0]

        ax_f.scatter(bins_mean_f_to_plot, hist_f_to_plot, alpha=0.9, s=10, c=species_color_map[species_name])
        # label=figure_utils.get_pretty_species_name(species_name)

        #print('mean', max(f_mean_mapgd))
        f_mean_mapgd_log10 = numpy.log10(f_mean_mapgd)
        f_mean_mapgd_log10_rescaled = (f_mean_mapgd_log10 - numpy.mean(f_mean_mapgd_log10)) / numpy.std(f_mean_mapgd_log10)
        hist_f_mean, bin_edges_f_mean = numpy.histogram(f_mean_mapgd_log10_rescaled, density=True, bins=20)
        bins_mean_f_mean = [0.5 * (bin_edges_f_mean[i] + bin_edges_f_mean[i+1]) for i in range(0, len(bin_edges_f_mean)-1 )]
        bins_mean_f_mean = numpy.asarray(bins_mean_f)

        hist_f_mean_to_plot = hist_f_mean[hist_f_mean>0]
        bins_mean_f_mean_to_plot = bins_mean_f_mean[hist_f_mean>0]

        ax_f_mean.scatter(bins_mean_f_mean_to_plot, hist_f_mean_to_plot, alpha=0.9, s=10, c=species_color_map[species_name])

        observed_prevalence_mapgd_no_zeros = observed_prevalence_mapgd[observed_prevalence_mapgd>0]
        observed_prevalence_mapgd_log10 = numpy.log10(observed_prevalence_mapgd_no_zeros)
        observed_prevalence_mapgd_log10_rescaled = (observed_prevalence_mapgd_log10 - numpy.mean(observed_prevalence_mapgd_log10)) / numpy.std(observed_prevalence_mapgd_log10)
        hist_observed_prevalence_mapgd, bin_edges_observed_prevalence_mapgd = numpy.histogram(observed_prevalence_mapgd_log10_rescaled, density=True, bins=20)
        bins_mean_observed_prevalence_mapgd = [0.5 * (bin_edges_observed_prevalence_mapgd[i] + bin_edges_observed_prevalence_mapgd[i+1]) for i in range(0, len(bin_edges_observed_prevalence_mapgd)-1 )]
        #ax_f_prevalence.scatter(bins_mean_observed_prevalence_mapgd, hist_observed_prevalence_mapgd, alpha=0.9, s=10, c=species_color_map[species_name])


        # taylors law

        #f_mean_mapgd_to_plot_test = f_mean_mapgd_to_plot[f_mean_mapgd_to_plot < 0.35]
        #f_var_mapgd_to_plot_test = f_var_mapgd_to_plot[f_mean_mapgd_to_plot < 0.35]

        ax_f_mean_vs_var.scatter(f_mean_mapgd_to_plot, f_var_mapgd_to_plot, alpha=1, s=12, c=species_color_map[species_name])
        means_all.extend(f_mean_mapgd_to_plot.tolist())
        variances_all.extend(f_var_mapgd_to_plot.tolist())

        set_n_non_zero_f = list(set(n_non_zero_f))
        set_n_non_zero_f.sort()

        max_n_to_plot = []
        prob_prevalence_to_plot = []
        #survival_array = [sum(n_non_zero_f>=i)/len(n_non_zero_f) for i in max_n_occurances_range]
        for n_i in max_n_occurances_range:

            if sum(n_non_zero_f==n_i) == 0:
                continue

            max_n_to_plot.append(n_i)
            prob_prevalence_to_plot.append(sum(n_non_zero_f==n_i)/len(n_non_zero_f))

        ax_f_prevalence.plot(max_n_to_plot, prob_prevalence_to_plot, lw=2, ls='-', marker='o', c=species_color_map[species_name], alpha=0.9)
        #ax_f_prevalence.plot(max_n_occurances_range, survival_array, lw=2, ls='-', alpha=0.8, c=species_color_map[species_name])

        #prevalence_range = numpy.logspace(-3, 0, num=1000, endpoint=True)
        #survival_array = [sum(observed_prevalence_mapgd>=i)/len(observed_prevalence_mapgd) for i in prevalence_range]
        #survival_array = numpy.asarray(survival_array)
        #ax_f_prevalence.plot(prevalence_range, survival_array, lw=2, ls='-', c=species_color_map[species_name], alpha=0.5)


        species_list_to_plot_good.append(species_name)

    ax_f_mean.set_xlabel('Rescaled ' + r'$\mathrm{log}_{10}$' + ' mean within-host\nallele frequency across hosts', fontsize=11)
    ax_f_mean.set_ylabel('Probability density', fontsize=12)
    #ax_f_mean.set_ylim([-0.02, 0.9])
    ax_f_mean.set_yscale('log', basey=10)
    ax_f_mean.xaxis.set_tick_params(labelsize=8)
    ax_f_mean.yaxis.set_tick_params(labelsize=8)

    if len(means_all) > 0:

        means_all = numpy.asarray(means_all)
        variances_all = numpy.asarray(variances_all)

        means_log10_all = numpy.log10(means_all)
        variances_log10_all = numpy.log10(variances_all)

        means_log10_all_test = means_log10_all[(means_log10_all > numpy.log10(4*(10**-2))) & (means_log10_all <= numpy.log10(f_mean_upper_cutoff_taylor)) ]
        variances_log10_all_test = variances_log10_all[(means_log10_all > numpy.log10(4*(10**-2))) & (means_log10_all <= numpy.log10(f_mean_upper_cutoff_taylor))]

        if len(means_log10_all_test) > 0:

            slope, intercept, r_value, p_value, std_err = stats.linregress(means_log10_all_test, variances_log10_all_test)

            lower_ci, upper_ci = prevalence_utils.calculate_slope_ci(means_log10_all_test, variances_log10_all_test, n=int(len(means_log10_all_test)))
            print('95% CI: ', lower_ci, upper_ci)

            # gamma AFD
            x_range = numpy.linspace(-4, 3, 10000)
            #k = 2.0
            k = slope
            k_digamma = special.digamma(k)
            k_trigamma = special.polygamma(1,k)
            gammalog = k*k_trigamma*x_range - numpy.exp(numpy.sqrt(k_trigamma)*x_range + k_digamma) - numpy.log(special.gamma(k)) + k*k_digamma + numpy.log10(numpy.exp(1))
            ax_f.plot(x_range, 10**gammalog, 'k', label='Gamma', lw=2)

            #x_log10_range =  numpy.linspace(min(means_log10_all) , max(means_log10_all) , 10000)
            x_log10_range =  numpy.linspace(min(means_log10_all_test) , max(means_log10_all_test) , 10000)
            y_log10_fit_range = (slope*x_log10_range + intercept)
            ax_f_mean_vs_var.plot(10**x_log10_range, 10**y_log10_fit_range, c='k', lw=2.5, linestyle='--', zorder=2, label=r'$\sigma^{{2}}_{{f}} \sim \bar{{f}}\,^{{{}}}$'.format(round(slope, 2)))


    ax_f.set_xlabel('Rescaled ' + r'$\mathrm{log}_{10}$' + ' within-host allele frequency', fontsize=11)
    ax_f.set_ylabel('Probability density', fontsize=12)
    #ax_f.set_ylim([-0.02, 1.23])
    ax_f.set_ylim([0.007, 1.04])
    ax_f.set_yscale('log', basey=10)
    ax_f.legend(loc='upper left', fontsize=11)
    ax_f.xaxis.set_tick_params(labelsize=8)
    ax_f.yaxis.set_tick_params(labelsize=8)

    # plot inequqlity
    x_inequality = numpy.logspace(numpy.log10(0.03), 0, num=1000, endpoint=True)
    y_inequality = (1-x_inequality)*x_inequality
    ax_f_mean_vs_var.plot(x_inequality, y_inequality, c='k', lw=2.5, linestyle=':', zorder=1, label='Max. ' + r'$\sigma_{f}^{2}$')

    ax_f_mean_vs_var.set_xscale('log', basex=10)
    ax_f_mean_vs_var.set_yscale('log', basey=10)
    ax_f_mean_vs_var.set_xlabel('Mean within-host allele\nfrequency across hosts, ' + r'$\bar{f}$', fontsize=12)
    ax_f_mean_vs_var.set_ylabel('Variance of within-host allele\nfrequencies across hosts, ' + r'$\sigma^{2}_{f}$', fontsize=11)
    ax_f_mean_vs_var.xaxis.set_tick_params(labelsize=8)
    ax_f_mean_vs_var.yaxis.set_tick_params(labelsize=7.5)
    ax_f_mean_vs_var.legend(loc='upper left', fontsize=11)


    ax_f_prevalence.set_xlabel('Number of hosts where allele is present (' + r'$f>0$' + ')', fontsize=11)
    ax_f_prevalence.set_ylabel('Fraction of sites', fontsize=12)
    ax_f_prevalence.xaxis.set_tick_params(labelsize=8)
    ax_f_prevalence.yaxis.set_tick_params(labelsize=8)
    #ax_f_prevalence.set_xscale('log', basex=10)
    ax_f_prevalence.set_yscale('log', basey=10)
    #ax_f_prevalence.set_ylim([-0.02, 1.02])


    fig.tight_layout()
    fig.subplots_adjust(wspace=0.34, hspace=0.28)
    fig.savefig("%sdiversity_summary_%s.pdf" % (config.analysis_directory, variant_type), format='pdf', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()



if __name__=='__main__':

    for variant_type in ['4D', '1D']:

        make_plot(variant_type)
