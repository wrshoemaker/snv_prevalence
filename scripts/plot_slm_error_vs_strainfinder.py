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
from collections import Counter

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import plot_utils
import prevalence_utils

from scipy.stats import gamma, gaussian_kde
import scipy.stats as stats

import calculate_predicted_prevalence_mapgd

iter = 100000
numpy.random.seed(123456789)
random.seed(123456789)
min_n_species = 5

#species_name = 'Bacteroides_xylanisolvens_57185'
clade_type = 'all'
variant_type = '4D'
pi_type = 'pi_include_boundary'
best = True

if best == False:
    best_status = '_best'
else:
    best_status = ''


species_color_map, ordered_species_list = plot_utils.get_species_color_map()




def get_strain_abundances(species_name, samples):

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    samples_to_keep = []
    richness_to_keep = []
    for sample in samples:

        intermediate_strain_filename = intermediate_strain_filename_template % sample

        if os.path.isfile(intermediate_strain_filename) == False:
            continue

        with open(intermediate_strain_filename, 'rb') as handle:
            b = pickle.load(handle)

        if species_name in b:

            abundances = b[species_name]
            abundances = numpy.asarray(abundances)

            samples_to_keep.append(sample)
            richness_to_keep.append(len(abundances))

    return samples_to_keep, richness_to_keep


def get_samples_in_analysis(species_name, clade_type, variant_type):

    frequency_dict = calculate_predicted_prevalence_mapgd.load_frequency_dict(species_name, clade_type)
    samples = []
    for key, value in frequency_dict[variant_type].items():
        samples.extend(value['samples'])

    return list(set(samples))


def make_strain_biodiversity_dict(clade_type, variant_type):

    strain_biodiversity_dict = {}

    for species_name in prevalence_utils.species_to_run:

        samples = get_samples_in_analysis(species_name, clade_type, variant_type)

        samples_strain, richnes_strains = get_strain_abundances(species_name, samples)
        richnes_strains = numpy.asarray(richnes_strains)

        fraction_samples_strain = sum(richnes_strains>1)/len(richnes_strains)
        mean_richness = numpy.mean(richnes_strains)

        strain_biodiversity_dict[species_name] = {}
        strain_biodiversity_dict[species_name]['fraction_samples_strain'] = fraction_samples_strain
        strain_biodiversity_dict[species_name]['mean_richness'] = mean_richness

    intermediate_filename_template = config.data_directory+"strain_biodiversity_%s_%s.dat"
    intermediate_filename = intermediate_filename_template % (clade_type, variant_type)

    sys.stderr.write("Saving strain biodiversity dict...\n")
    with open(intermediate_filename, 'wb') as handle:
        pickle.dump(strain_biodiversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_strain_biodiversity_dict(clade_type, variant_type):

    intermediate_filename_template = config.data_directory+"strain_biodiversity_%s_%s.dat"
    intermediate_filename = intermediate_filename_template % (clade_type, variant_type)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)

    return b




def make_error_dict():

    sys.stderr.write("Generating error dict...\n")

    prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()
    species_list = list(prevalence_dict_mapgd.keys())
    species_list.sort()

    error_n_dict = {}

    for species_name in species_list:

        predicted_prevalence = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm%s' % best_status]
        predicted_prevalence = numpy.asarray(predicted_prevalence)

        observed_prevalence = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm%s' % best_status]
        observed_prevalence = numpy.asarray(observed_prevalence)

        n_non_zero_frequencies = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
        n_non_zero_frequencies = numpy.asarray(n_non_zero_frequencies)

        predicted_prevalence_no_zeros = predicted_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]
        observed_prevalence_no_zeros = observed_prevalence[(observed_prevalence>0) & (predicted_prevalence>0) ]
        n_non_zero_frequencies_no_zeros = n_non_zero_frequencies[(observed_prevalence>0) & (predicted_prevalence>0)]

        if len(observed_prevalence_no_zeros) == 0:
            continue

        n_non_zero_frequencies_no_zeros_set = list(set(n_non_zero_frequencies_no_zeros.tolist()))

        error = numpy.absolute(observed_prevalence_no_zeros - predicted_prevalence_no_zeros) / observed_prevalence_no_zeros
        #mean_error = numpy.mean(numpy.log10(error))
        mean_error = numpy.mean(error)

        pi_dict = calculate_predicted_prevalence_mapgd.load_pi_dict(species_name)
        samples = list(pi_dict[variant_type].keys())

        samples_strains, richnes_strains = get_strain_abundances(species_name, samples)

        richnes_strains = numpy.asarray(richnes_strains)

        fraction_samples_strain = sum(richnes_strains>1)/len(richnes_strains)
        mean_richness = numpy.mean(richnes_strains)

        error_n_dict[species_name] = {}

        error_n_dict[species_name]['mean_richness'] = mean_richness
        error_n_dict[species_name]['fraction_samples_strain'] = fraction_samples_strain
        error_n_dict[species_name]['mean_error'] = mean_error
        error_n_dict[species_name]['prevalence_error'] = {}

        error_n_dict[species_name]['error'] = error.tolist()
        error_n_dict[species_name]['n_non_zero_frequencies'] = n_non_zero_frequencies.tolist()
        error_n_dict[species_name]['observed_prevalence_no_zeros'] = observed_prevalence_no_zeros.tolist()
        error_n_dict[species_name]['error'] = error.tolist()

        for n in n_non_zero_frequencies_no_zeros_set:

            if sum(n_non_zero_frequencies_no_zeros==n) < 3:
                continue

            #mean_error_n = numpy.mean(numpy.log10(error[n_non_zero_frequencies_no_zeros==n]))
            mean_error_n = numpy.mean(error[n_non_zero_frequencies_no_zeros==n])

            error_n_dict[species_name]['prevalence_error'][n] = mean_error_n


    intermediate_filename_template = config.data_directory+"error_%s_%s%s.dat"
    intermediate_filename = intermediate_filename_template % (clade_type, variant_type, best_status)

    sys.stderr.write("Saving strain vs. error dict...\n")
    with open(intermediate_filename, 'wb') as handle:
        pickle.dump(error_n_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_error_dict():

    intermediate_filename_template = config.data_directory+"error_%s_%s%s.dat"
    intermediate_filename = intermediate_filename_template % (clade_type, variant_type, best_status)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)

    return b



def make_permutation_dict():

    sys.stderr.write("Generating permutation dict...\n")
    error_n_dict = load_error_dict()
    species_list = list(error_n_dict.keys())
    min_prevalence_all = []
    max_prevalence_all = []
    max_strain_fraction = max([error_n_dict[species_name]['fraction_samples_strain'] for species_name in species_list])
    # get all number of occurances
    for species_name in species_list:
        min_prevalence_all.append(min(error_n_dict[species_name]['observed_prevalence_no_zeros']))
        max_prevalence_all.append(max(error_n_dict[species_name]['observed_prevalence_no_zeros']))

    min_prevalence = max(min_prevalence_all)
    #max_prevalence = min(max_prevalence_all)
    max_prevalence = max_strain_fraction

    prevalence_species_cutoff_dict = {}
    for species_name in species_list:
        observed_prevalence_no_zeros = numpy.asarray(error_n_dict[species_name]['observed_prevalence_no_zeros'])
        error = numpy.asarray(error_n_dict[species_name]['error'])

        prevalence_species_cutoff_dict[species_name] = {}
        observed_prevalence_no_zeros_to_keep = observed_prevalence_no_zeros[(observed_prevalence_no_zeros>=min_prevalence) & (observed_prevalence_no_zeros <= max_prevalence)]
        error_to_keep = error[(observed_prevalence_no_zeros>=min_prevalence) & (observed_prevalence_no_zeros <= max_prevalence)]

        prevalence_species_cutoff_dict[species_name]['observed_prevalence_no_zeros'] = observed_prevalence_no_zeros_to_keep
        prevalence_species_cutoff_dict[species_name]['error'] = error_to_keep


    # Examine the range of prevalence values that all species have


    min_prevalence_range = numpy.logspace(numpy.log10(min_prevalence), numpy.log10(max_prevalence), num=100, endpoint=False, base=10.0)

    # make null distribution
    error_prevalence_dict_species = {}
    for species_name in species_list:
        error_prevalence_dict_species[species_name] = []

    rho_dict = {}
    species_samples_prevalence_dict = {}
    for min_prevalence_ in min_prevalence_range:

        mean_error_all = []
        strain_fraction_all = []
        species_all = []

        for species_name in species_list:

            observed_prevalence  = prevalence_species_cutoff_dict[species_name]['observed_prevalence_no_zeros']
            error = prevalence_species_cutoff_dict[species_name]['error']

            observed_prevalence_to_keep = observed_prevalence[(observed_prevalence>=min_prevalence_) & (observed_prevalence <= max_prevalence)]
            error_to_keep = error[(observed_prevalence>=min_prevalence_) & (observed_prevalence <= max_prevalence)]

            if len(error_to_keep) < 10:
                continue

            mean_error_all.append(numpy.mean(error_to_keep))
            strain_fraction_all.append(error_n_dict[species_name]['fraction_samples_strain'])
            species_all.append(species_name)

        # at least three observations
        if len(strain_fraction_all) < min_n_species:
            continue

        # observed
        rho = numpy.corrcoef(numpy.log10(mean_error_all), strain_fraction_all)[0,1]
        rho_dict[min_prevalence_] = rho
        species_samples_prevalence_dict[min_prevalence_] = species_all

        # add data to dict for permutation

        for species_name in species_all:
            error_prevalence_dict_species[species_name].append(min_prevalence_)


    null_error_prevalence_dict_species = {}
    for min_prevalence_ in min_prevalence_range:
        null_error_prevalence_dict_species[min_prevalence_] = {}

    for i in range(iter):

        for species_name in species_list:

            observed_prevalence  = prevalence_species_cutoff_dict[species_name]['observed_prevalence_no_zeros']
            error = prevalence_species_cutoff_dict[species_name]['error']
            error_permute = numpy.random.permutation(error)

            for min_prevalence_ in min_prevalence_range:

                if min_prevalence_ not in error_prevalence_dict_species[species_name]:
                    continue

                observed_prevalence_to_keep = observed_prevalence[(observed_prevalence>=min_prevalence_) & (observed_prevalence <= max_prevalence)]
                error_permute_to_keep = error_permute[(observed_prevalence>=min_prevalence_) & (observed_prevalence <= max_prevalence)]

                if species_name not in null_error_prevalence_dict_species[min_prevalence_]:
                    null_error_prevalence_dict_species[min_prevalence_][species_name] = []
                null_error_prevalence_dict_species[min_prevalence_][species_name].append(numpy.mean(error_permute_to_keep))

    min_prevalence_to_keep = list(null_error_prevalence_dict_species.keys())
    rho_null_dict = {}
    for min_prevalence_ in min_prevalence_to_keep:

        species_min_prevalence = list(null_error_prevalence_dict_species[min_prevalence_].keys())

        if len(species_min_prevalence) < min_n_species:
            continue

        rho_null_dict[min_prevalence_] = []
        for i_ in range(iter):
            null_strain_fraction_min_prevalence_all = [error_n_dict[s]['fraction_samples_strain'] for s in species_min_prevalence]
            null_error_min_prevalence_all = [null_error_prevalence_dict_species[min_prevalence_][s][i_] for s in species_min_prevalence]
            rho_null = numpy.corrcoef(numpy.log10(null_error_min_prevalence_all), null_strain_fraction_min_prevalence_all)[0,1]
            rho_null_dict[min_prevalence_].append(rho_null)


    permutation_dict = {}

    min_prevalence_to_save = list(rho_null_dict.keys())
    min_prevalence_to_save.sort()
    for min_prevalence_i in min_prevalence_to_save:

        rho_null_dict_n_i = numpy.asarray(rho_null_dict[min_prevalence_i])
        rho_null_dict_n_i = numpy.sort(rho_null_dict_n_i)
        lower_ci = rho_null_dict_n_i[int(0.025*iter)]
        upper_ci = rho_null_dict_n_i[int(0.975*iter)]

        permutation_dict[min_prevalence_i] = {}
        permutation_dict[min_prevalence_i]['observed'] = rho_dict[min_prevalence_i]
        permutation_dict[min_prevalence_i]['lower_ci'] = lower_ci
        permutation_dict[min_prevalence_i]['upper_ci'] = upper_ci
        permutation_dict[min_prevalence_i]['species'] = species_samples_prevalence_dict[min_prevalence_i]

        print(min_prevalence_i, len(species_samples_prevalence_dict[min_prevalence_i]))


    sys.stderr.write("Saving permutation dict...\n")
    intermediate_filename = "%sslm_error_strain_corr_%s%s.dat" % (config.data_directory, clade_type, best_status)
    with open(intermediate_filename, 'wb') as handle:
        pickle.dump(permutation_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_permutation_dict():

    sys.stderr.write("Loading permutation dict...\n")
    intermediate_filename = "%sslm_error_strain_corr_%s%s.dat" % (config.data_directory, clade_type, best_status)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b



#make_error_dict()
#make_permutation_dict()

error_n_dict = load_error_dict()
species_name_all = list(error_n_dict.keys())

#slope, intercept, r_value, p_value, std_err = stats.linregress(fraction_samples_strain_all, all_error_all)

permutation_dict = load_permutation_dict()
n_to_plot = list(permutation_dict.keys())
n_to_plot.sort()

lower_ci = [permutation_dict[n_]['lower_ci'] for n_ in n_to_plot]
upper_ci = [permutation_dict[n_]['upper_ci'] for n_ in n_to_plot]
observed = [permutation_dict[n_]['observed'] for n_ in n_to_plot]



# get examples for regression
prevalence_error_31 = []
prevalence_error_1 = []
fraction_samples_strain_31 = []
fraction_samples_strain_1 = []
species_31 = []
species_1 = []
for species_name in error_n_dict.keys():

    if 1 in error_n_dict[species_name]['prevalence_error']:
        prevalence_error_1.append(error_n_dict[species_name]['prevalence_error'][1])
        fraction_samples_strain_1.append(error_n_dict[species_name]['fraction_samples_strain'])
        species_1.append(species_name)


    if 35 in error_n_dict[species_name]['prevalence_error']:
        prevalence_error_31.append(error_n_dict[species_name]['prevalence_error'][35])
        fraction_samples_strain_31.append(error_n_dict[species_name]['fraction_samples_strain'])
        species_31.append(species_name)





#fig, ax = plt.subplots(figsize=(4,4))
gs = gridspec.GridSpec(nrows=1, ncols=3)
fig = plt.figure(figsize = (14, 4))
ax_strain = fig.add_subplot(gs[0, 0])
ax_ex = fig.add_subplot(gs[0, 1])
ax_corr =fig.add_subplot(gs[0, 2])

ax_strain.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_strain.transAxes)
ax_ex.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_ex.transAxes)
ax_corr.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_corr.transAxes)


# plot strain bar plot
species_fraction_strains_tuple = [(s, error_n_dict[s]['fraction_samples_strain']) for s in error_n_dict.keys()]
species_fraction_strains_tuple.sort(key=lambda tup: tup[1])
species_fraction_strains_tuple = species_fraction_strains_tuple[::-1]

species_to_plot = [s[0] for s in species_fraction_strains_tuple]
strain_fraction_to_plot = [s[1] for s in species_fraction_strains_tuple]

species_to_plot_pretty = [figure_utils.get_pretty_species_name(s) for s in species_to_plot]
colors_hosts = [species_color_map[s] for s in species_to_plot]
ax_strain.barh(species_to_plot_pretty, strain_fraction_to_plot, height=0.8, align='center', color=colors_hosts)
ax_strain.set_xlabel('Strain prevalence\n(fraction of hosts with strain structure)', fontsize=11)
ax_strain.xaxis.set_tick_params(labelsize=8)
ax_strain.yaxis.set_tick_params(labelsize=9)
ax_strain.set_ylim([-0.6, len(species_to_plot)-0.3])


#sorted_samples = [s[1] for s in sorted(zip(good_species_list, number_samples), key = lambda t: t[1])][::-1]

#good_species_list_sorted_samples_pretty = [figure_utils.get_pretty_species_name(s) for s in good_species_list_sorted_samples][::-1]

#error_n_dict[s]['fraction_samples_strain']




ax_ex.scatter(-100, -100, color='k', linewidth=2, facecolors='white', s=60, label='SNVs present in one host')
ax_ex.scatter(-100, -100, color='k', s=60, label='SNVs present in 35 hosts')
ax_ex.set_xlim([0.05, 0.7])
ax_ex.set_ylim([0.08, 0.92])
ax_ex.set_xlabel('Strain prevalence', fontsize=11)
ax_ex.set_ylabel('Mean relative error of SLM', fontsize=11)
ax_ex.legend(loc="upper right", fontsize=7)




slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(fraction_samples_strain_1, prevalence_error_1)
slope_31, intercept_31, r_value_31, p_value_31, std_err_31 = stats.linregress(fraction_samples_strain_31, prevalence_error_31)

x_range =  numpy.linspace((0.1, 0.6) , 10000)


y_fit_1 = (slope_1*x_range + intercept_1)
y_fit_31 = (slope_31*x_range + intercept_31)

ax_ex.plot(x_range, y_fit_1, c='k', lw=2.5, linestyle='--', zorder=2)
ax_ex.plot(x_range, y_fit_31, c='k', lw=2.5, linestyle='-', zorder=2)
ax_ex.xaxis.set_tick_params(labelsize=8)
ax_ex.yaxis.set_tick_params(labelsize=8)
ax_ex.text(0.83,  0.25, r'$\rho_{1} =$' + str( round(r_value_1, 3) ), fontsize=8, color='k', ha='center', va='center', transform=ax_ex.transAxes  )
ax_ex.text(0.578 , 0.25, r'$P_{1} =$' + str( round(p_value_1, 3) ), fontsize=8, color='k', ha='center', va='center', transform=ax_ex.transAxes  )
ax_ex.text(0.85, 0.2, r'$\rho_{35} =$' + str(round(r_value_31, 3)), fontsize=8, color='k', ha='center', va='center', transform=ax_ex.transAxes  )
ax_ex.text(0.6, 0.2, r'$P_{35}=$' + str(round(p_value_31, 4) ), fontsize=8, color='k', ha='center', va='center', transform=ax_ex.transAxes  )

ax_ex.set_ylim([0.0, 0.95])

#ax_ex.text(0.6,0.2, 'ddf', fontsize=11, color='k', ha='center', va='center', transform=ax_ex.transAxes  )


ax_corr.fill_between(n_to_plot, lower_ci, upper_ci, alpha=0.6, color='grey', zorder=1, label='95% CI')
ax_corr.plot(n_to_plot, lower_ci, color ='k', lw=1.5, zorder=1)
ax_corr.plot(n_to_plot, upper_ci, color ='k', lw=1.5, zorder=1 )
ax_corr.plot(n_to_plot, observed, marker="o", color='k', lw=2, label="Observed")
ax_corr.set_xlim([min(n_to_plot), max(n_to_plot)])
ax_corr.legend(loc='upper left', fontsize=8)
ax_corr.set_xlabel('Observed SNV prevalence', fontsize=11)
ax_corr.set_ylabel('Correlation between\nSLM MRE and strain prevalence', fontsize=11)
ax_corr.axhline(0, lw=1.5, ls='--',color='k', zorder=1)
ax_corr.xaxis.set_tick_params(labelsize=8)
ax_corr.yaxis.set_tick_params(labelsize=6.5)

fig.tight_layout()
fig.subplots_adjust(hspace=0.1, wspace=0.28)
fig.savefig("%sslm_error_vs_strainfinder%s.png" % (config.analysis_directory, best_status), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
plt.close()
