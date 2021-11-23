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
import parse_HMP_data

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from scipy.stats import gamma, gaussian_kde
import scipy.stats as stats

import prevalence_utils

import calculate_predicted_occupancy

numpy.random.seed(123456789)

min_sample_size = config.between_host_min_sample_size
low_divergence_threshold = config.between_low_divergence_threshold

allowed_variant_types = set(['1D','4D'])

good_species_list = parse_midas_data.parse_good_species_list()
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])
clade_types = ['all','largest_clade']

clade_type_label_dict = {'all': 'All hosts', 'largest_clade': 'Largest clade'}


prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()



def generate_date():

    out_dict = {}

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    for species_name in good_species_list:

        # get largest clade sampels and all others

        #sys.stderr.write("Loading whitelisted genes...\n")
        #core_genes = core_gene_utils.parse_core_genes(species_name)

        #snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        #if os.path.isfile(snp_file_path) == False:
        #    continue

        #snp_file =  bz2.BZ2File(snp_file_path, "r")
        #line = snp_file.readline() # header
        #items = line.split()[1:]

        #samples = numpy.array([item.strip() for item in items])
        #samples_in_snp_and_pi_dict = numpy.intersect1d(samples, samples_from_dict)

        # get samples from largest clade
        snp_samples = diversity_utils.calculate_haploid_samples(species_name)
        if len(snp_samples) < min_sample_size:
            continue
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        # get clade idxs
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = numpy.array(dummy_samples)
        substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        coarse_grained_idxs, coarse_grained_cluster_list = clade_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold) # NRG: what is this returning?

        coarse_grained_samples = snp_samples[coarse_grained_idxs]
        clade_sets = clade_utils.load_manual_clades(species_name)

        if len(clade_sets)==0:
            continue

        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)
        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])
        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]
        largest_clade_set = set(largest_clade_samples)

        #snp_file.close()


        #####


        snp_file_path = "%ssnps/%s/snps_summary.txt" % (config.data_directory, species_name)

        if os.path.isfile(snp_file_path) == False:
            continue

        samples = []

        for line in open(snp_file_path):

            line = line.strip().split('\t')
            sample = line[0]
            samples.append(sample)


        samples = numpy.array(samples)

        richness = []
        diversity = []
        evenness = []
        dominance = []

        richness_largest_clade = []
        diversity_largest_clade = []
        evenness_largest_clade = []
        dominance_largest_clade = []

        #n_samples_with_strains = 0
        #n_samples_with_strains_largest_clade = 0
        # get strain data
        for sample in samples:

            intermediate_strain_filename = intermediate_strain_filename_template % sample

            if os.path.isfile(intermediate_strain_filename) == False:
                continue

            with open(intermediate_strain_filename, 'rb') as handle:
                b = pickle.load(handle)

            if species_name in b:

                abundances = b[species_name]
                abundances = numpy.asarray(abundances)

                richness.append(len(abundances))

                if sample in largest_clade_set:
                    richness_largest_clade.append(len(abundances))

                if len(abundances) > 1:

                    #n_samples_with_strains += 1
                    H = -1*sum(abundances*numpy.log(abundances))
                    J = H/numpy.log(len(abundances))

                    diversity.append(H)
                    evenness.append(J)
                    dominance.append(max(abundances))

                    if sample in largest_clade_set:

                        diversity_largest_clade.append(diversity)
                        evenness_largest_clade.append(evenness)
                        dominance_largest_clade.append(max(abundances))




        richness = numpy.asarray(richness)
        richness_largest_clade = numpy.asarray(richness_largest_clade)

        # skip if less than three observations
        if len(richness) < 3:
            continue

        out_dict[species_name] = {}
        out_dict[species_name]['all'] = {}
        out_dict[species_name]['all']['fraction_samples_strain'] = sum(richness>1)/len(richness)


        if len(diversity)>=3:

            out_dict[species_name]['all']['mean_richness'] = numpy.mean(richness)
            out_dict[species_name]['all']['mean_diversity'] = numpy.mean(diversity)
            out_dict[species_name]['all']['mean_evenness'] = numpy.mean(evenness)
            out_dict[species_name]['all']['mean_dominance'] = numpy.mean(dominance)

        #if len(diversity_largest_clade) < 3:
        #    continue

        if len(richness_largest_clade) < 3:
            continue

        out_dict[species_name]['largest_clade'] = {}
        out_dict[species_name]['largest_clade']['fraction_samples_strain'] = sum(richness_largest_clade>1)/len(richness_largest_clade)

        if len(diversity_largest_clade)>=3:

            out_dict[species_name]['largest_clade']['mean_richness'] = numpy.mean(richness_largest_clade)
            out_dict[species_name]['largest_clade']['mean_diversity'] = numpy.mean(diversity_largest_clade)
            out_dict[species_name]['largest_clade']['mean_evenness'] = numpy.mean(evenness_largest_clade)
            out_dict[species_name]['largest_clade']['mean_dominance'] = numpy.mean(dominance_largest_clade)



    output_filename = config.data_directory+"strain_diversity.dat"

    with open(output_filename, 'wb') as handle:
        pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




output_filename = config.data_directory+"strain_diversity.dat"

with open(output_filename, 'rb') as handle:
    strain_diversity_dict = pickle.load(handle)


#fig, ax = plt.subplots(figsize=(4,4))

fig = plt.figure(figsize = (8, 4)) #

ax_mae = plt.subplot2grid((2, 3), (0, 0))
ax_me = plt.subplot2grid((2, 3), (1, 0))

ax_mae.text(-0.1, 1.07, prevalence_utils.sub_plot_labels[0], fontsize=8, fontweight='bold', ha='center', va='center', transform=ax_mae.transAxes)
ax_me.text(-0.1, 1.07, prevalence_utils.sub_plot_labels[1], fontsize=8, fontweight='bold', ha='center', va='center', transform=ax_me.transAxes)

ax_strain = plt.subplot2grid((2, 3), (0,1), colspan=1)
ax_richness = plt.subplot2grid((2, 3), (1,1), colspan=1)

ax_strain.text(-0.1, 1.07, prevalence_utils.sub_plot_labels[2], fontsize=8, fontweight='bold', ha='center', va='center', transform=ax_strain.transAxes)
ax_richness.text(-0.1, 1.07, prevalence_utils.sub_plot_labels[3], fontsize=8, fontweight='bold', ha='center', va='center', transform=ax_richness.transAxes)

ax_dominance = plt.subplot2grid((2, 3), (0,2), colspan=1)
ax_evenness = plt.subplot2grid((2, 3), (1,2), colspan=1)

ax_dominance.text(-0.1, 1.07, prevalence_utils.sub_plot_labels[4], fontsize=8, fontweight='bold', ha='center', va='center', transform=ax_dominance.transAxes)
ax_evenness.text(-0.1, 1.07, prevalence_utils.sub_plot_labels[5], fontsize=8, fontweight='bold', ha='center', va='center', transform=ax_evenness.transAxes)



species_list = list(prevalence_dict.keys())

mae_all = []
mae_largest_clade = []

me_all = []
me_largest_clade = []

for species_name in species_list:

    mae_list = []
    me_list = []

    for clade_type in clade_types:


        mae = prevalence_dict[species_name][clade_type]['4D']['MRE']
        me = prevalence_dict[species_name][clade_type]['4D']['ME']

        mae_list.append(mae)
        me_list.append(me)

    mae_all.append(mae_list[0])
    mae_largest_clade.append(mae_list[1])


    me_all.append(me_list[0])
    me_largest_clade.append(me_list[1])


    ax_mae.plot([0,1], mae_list, c='k', alpha=0.7, zorder=1)
    ax_me.plot([0,1], me_list, c='k', alpha=0.7, zorder=1)

    ax_mae.scatter([0,1], mae_list, alpha=0.7, s=40, c='k', zorder=2)
    ax_me.scatter([0,1], me_list, alpha=0.7, s=40, c='k', zorder=2)



t_mae, p_mae = prevalence_utils.permutational_two_sample_t_test_equal_n(mae_all, mae_largest_clade)

t_me, p_me = prevalence_utils.permutational_two_sample_t_test_equal_n(me_all, me_largest_clade)


ax_mae.text(0.18,0.93, r'$t = $' + str(round(t_mae, 2)), fontsize=8, color='k', ha='center', va='center', transform=ax_mae.transAxes  )
ax_mae.text(0.2,0.85, r'$P = $' + str(round(p_mae, 4)), fontsize=8, color='k', ha='center', va='center', transform=ax_mae.transAxes  )



ax_me.text(0.8,0.3, r'$t = $' + str(round(t_me, 2)), fontsize=8, color='k', ha='center', va='center', transform=ax_me.transAxes  )
ax_me.text(0.82,0.19, r'$P = $' + str(round(p_me, 4)), fontsize=8, color='k', ha='center', va='center', transform=ax_me.transAxes  )



ax_mae.set_xticks([0,1])
ax_me.set_xticks([0,1])

ax_mae.set_xlim([-0.2,1.2])
ax_me.set_xlim([-0.2,1.2])

ax_mae.set_xticklabels(['All samples', 'Largest clade'], fontsize=8)
ax_me.set_xticklabels(['All samples', 'Largest clade'], fontsize=8)

ax_mae.set_ylabel('Mean relative error', fontsize=8)
ax_me.set_ylabel('Mean error', fontsize=8)

#ax_mae.set_yscale('log', basey=10)

ax_mae.yaxis.set_tick_params(labelsize=6)
ax_me.yaxis.set_tick_params(labelsize=6)



legend_elements = [Line2D([0], [0], marker='o', color='w', label='Species',
                          markerfacecolor='k', markersize=8)]

#ax_mae.legend(handles=legend_elements, loc='upper left', fontsize='small')




indices_to_plot = {}

for index in ['richness', 'evenness', 'dominance', 'diversity', 'fraction_samples_strain']:
    indices_to_plot[index] = {}
    indices_to_plot[index]['mae'] = []
    indices_to_plot[index]['index_values'] = []


for species_name in good_species_list:

    if (species_name not in prevalence_dict) or (species_name not in strain_diversity_dict):
        continue


    mae = prevalence_dict[species_name]['all']['4D']['MRE']


    if 'fraction_samples_strain' in strain_diversity_dict[species_name]['all']:
        fraction_samples_strain = strain_diversity_dict[species_name]['all']['fraction_samples_strain']
        indices_to_plot['fraction_samples_strain']['mae'].append(mae)
        indices_to_plot['fraction_samples_strain']['index_values'].append(fraction_samples_strain)


    if 'mean_richness' in strain_diversity_dict[species_name]['all']:
        richness = strain_diversity_dict[species_name]['all']['mean_richness']
        indices_to_plot['richness']['mae'].append(mae)
        indices_to_plot['richness']['index_values'].append(richness)

    if 'mean_diversity' in strain_diversity_dict[species_name]['all']:
        diversity = strain_diversity_dict[species_name]['all']['mean_diversity']
        indices_to_plot['diversity']['mae'].append(mae)
        indices_to_plot['diversity']['index_values'].append(diversity)


    if 'mean_evenness' in strain_diversity_dict[species_name]['all']:
        evenness = strain_diversity_dict[species_name]['all']['mean_evenness']
        indices_to_plot['evenness']['mae'].append(mae)
        indices_to_plot['evenness']['index_values'].append(evenness)


    if 'mean_dominance' in strain_diversity_dict[species_name]['all']:
        dominance = strain_diversity_dict[species_name]['all']['mean_dominance']
        indices_to_plot['dominance']['mae'].append(mae)
        indices_to_plot['dominance']['index_values'].append(dominance)



ax_strain.scatter(indices_to_plot['fraction_samples_strain']['mae'], indices_to_plot['fraction_samples_strain']['index_values'], c='k', alpha=0.7, s=60, label='Species')
ax_strain.set_xlabel('Mean relative error of the gamma', fontsize=8)
ax_strain.set_ylabel('Fraction hosts with > 1 strain', fontsize=8)
ax_strain.set_xlim([min(indices_to_plot['fraction_samples_strain']['mae'])*0.8, max(indices_to_plot['fraction_samples_strain']['mae'])*1.1])
ax_strain.set_ylim([min(indices_to_plot['fraction_samples_strain']['index_values'])*0.8, max(indices_to_plot['fraction_samples_strain']['index_values'])*1.1])
ax_strain.xaxis.set_tick_params(labelsize=6)
ax_strain.yaxis.set_tick_params(labelsize=6)

slope, intercept, r_value, p_value, std_err = stats.linregress(indices_to_plot['fraction_samples_strain']['mae'], indices_to_plot['fraction_samples_strain']['index_values'])

if p_value >= 0.05:
    ax_strain.text(0.85,0.9, r'$P \nless 0.05$', fontsize=8, color='k', ha='center', va='center', transform=ax_strain.transAxes  )





ax_richness.scatter(indices_to_plot['richness']['mae'], indices_to_plot['richness']['index_values'], c='k', alpha=0.7, s=60, label='One species')
ax_richness.set_xlabel('Mean relative error of the gamma', fontsize=8)
ax_richness.set_ylabel('Mean strain richness of\nhosts with $>1$ strain', fontsize=8)
ax_richness.set_xlim([min(indices_to_plot['richness']['mae'])*0.8, max(indices_to_plot['richness']['mae'])*1.1])
ax_richness.set_ylim([min(indices_to_plot['richness']['index_values'])*0.8, max(indices_to_plot['richness']['index_values'])*1.1])
ax_richness.xaxis.set_tick_params(labelsize=6)
ax_richness.yaxis.set_tick_params(labelsize=6)

slope, intercept, r_value, p_value, std_err = stats.linregress(indices_to_plot['richness']['mae'], indices_to_plot['richness']['index_values'])

if p_value >= 0.05:
    ax_richness.text(0.85,0.9, r'$P \nless 0.05$', fontsize=8, color='k', ha='center', va='center', transform=ax_richness.transAxes  )



ax_dominance.scatter(indices_to_plot['dominance']['mae'], indices_to_plot['dominance']['index_values'], c='k', alpha=0.7, s=60, label='One species')
ax_dominance.set_xlabel('Mean relative error of the gamma', fontsize=8)
ax_dominance.set_ylabel('Mean strain dominance of\nhosts with $>1$ strain', fontsize=8)
ax_dominance.axhline(y=0.5, color='k', linestyle=':', lw = 3, zorder=1, label="Minimum dominance")
ax_dominance.set_xlim([min(indices_to_plot['dominance']['mae'])*0.8, max(indices_to_plot['dominance']['mae'])*1.1])
ax_dominance.set_ylim([0.48, max(indices_to_plot['dominance']['index_values'])*1.1])
ax_dominance.legend(loc="upper right", prop={'size': 6})
ax_dominance.xaxis.set_tick_params(labelsize=6)
ax_dominance.yaxis.set_tick_params(labelsize=6)
slope, intercept, r_value, p_value, std_err = stats.linregress(indices_to_plot['dominance']['mae'], indices_to_plot['dominance']['index_values'])

if p_value >= 0.05:
    ax_dominance.text(0.14,0.9, r'$P \nless 0.05$', fontsize=8, color='k', ha='center', va='center', transform=ax_dominance.transAxes  )



ax_evenness.scatter(indices_to_plot['evenness']['mae'], indices_to_plot['evenness']['index_values'], c='k', alpha=0.7, s=60, label='Species')
ax_evenness.set_xlabel('Mean relative error of the gamma', fontsize=8)
ax_evenness.set_ylabel('Mean strain evenness of\nhosts with $>1$ strain', fontsize=8)
ax_evenness.set_xlim([min(indices_to_plot['evenness']['mae'])*0.8, max(indices_to_plot['evenness']['mae'])*1.1])
ax_evenness.set_ylim([min(indices_to_plot['evenness']['index_values'])*0.9, max(indices_to_plot['evenness']['index_values'])*1.1])
ax_evenness.xaxis.set_tick_params(labelsize=6)
ax_evenness.yaxis.set_tick_params(labelsize=6)
slope, intercept, r_value, p_value, std_err = stats.linregress(indices_to_plot['evenness']['mae'], indices_to_plot['evenness']['index_values'])

if p_value >= 0.05:
    ax_evenness.text(0.85,0.9, r'$P \nless 0.05$', fontsize=8, color='k', ha='center', va='center', transform=ax_evenness.transAxes  )



fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%sprevalence_mse_strain.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
