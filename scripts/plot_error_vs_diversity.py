from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path
import itertools


import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
import scipy.stats as stats
import scipy.special as special

#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.gridspec as gridspec
#from matplotlib.lines import Line2D

import calculate_predicted_prevalence_mapgd

#import plot_utils


error_vs_diversity_dict_template = config.data_directory+"error_vs_diversity_dict/%s.dat"



prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all(test=True)
species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()

variant_types = ['4D']
clade_type = 'all'
pi_type = 'pi_include_boundary'
min_n_sites = 100

iter = 100

def get_frequency_dict_mapgd(species_name, clade_type='all'):

    intermediate_filename_template = config.data_directory+"frequency_dict_mapgd/%s_%s.dat"
    intermediate_filename = intermediate_filename_template % (species_name, clade_type)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b




def make_error_vs_diversity_dict(species_name):

    sys.stderr.write("Running %s...\n" % species_name)

    error_vs_diversity_dict = {}

    frequency_dict_mapgd = get_frequency_dict_mapgd(species_name)

    for variant_type in variant_types:

        sys.stderr.write("Running %s sites...\n" % variant_type)

        error_vs_diversity_dict[variant_type] = {}

        frequency_dict_mapgd_keys = list(frequency_dict_mapgd[variant_type].keys())
        # rename keys to get rid of "L"
        for k in frequency_dict_mapgd_keys:
            #k_new = (k[0], int(k[1]))
            k_new = "%s_%s" % (k[0], str(int(k[1])))
            frequency_dict_mapgd[variant_type][k_new] = frequency_dict_mapgd[variant_type].pop(k)

        observed_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['observed_prevalence_mapgd_slm']
        observed_prevalence_mapgd = numpy.asarray(observed_prevalence_mapgd)

        predicted_prevalence_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd_slm']
        predicted_prevalence_mapgd = numpy.asarray(predicted_prevalence_mapgd)

        n_non_zero_f = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['n_non_zero_frequencies']
        n_non_zero_f = numpy.asarray(n_non_zero_f)

        f_max_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_max_mapgd']
        f_max_mapgd = numpy.asarray(f_max_mapgd)

        sites_mapgd = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['sites']
        sites_mapgd = numpy.asarray(sites_mapgd)

        idx_ = (observed_prevalence_mapgd>0) & (predicted_prevalence_mapgd>0) & (f_max_mapgd<1)
        observed_prevalence_mapgd = observed_prevalence_mapgd[idx_]
        predicted_prevalence_mapgd = predicted_prevalence_mapgd[idx_]
        n_non_zero_f = n_non_zero_f[idx_]
        sites_mapgd = sites_mapgd[idx_]

        error = numpy.absolute(observed_prevalence_mapgd-predicted_prevalence_mapgd)/observed_prevalence_mapgd

        site_mapgd_to_keep = []
        samples_all = []
        error_to_keep = []
        for s_idx, s in enumerate(sites_mapgd):
            s_tuple = "%s_%s" % (s[0], str(int(s[1])))
            if s_tuple in frequency_dict_mapgd[variant_type]:
                site_mapgd_to_keep.append(s_tuple)
                samples_all.extend(frequency_dict_mapgd[variant_type][s_tuple]['samples'])
                error_to_keep.append(error[s_idx])

        site_mapgd_to_keep = numpy.asarray(site_mapgd_to_keep)
        error_to_keep = numpy.asarray(error_to_keep)
        #samples_all = list(set(samples_all))
        frequency_sample_dict = {}
        for s in site_mapgd_to_keep:

            samples_s = frequency_dict_mapgd[variant_type][s]['samples']
            samples_non_zero_s = frequency_dict_mapgd[variant_type][s]['samples_non_zero']
            frequencies_s = frequency_dict_mapgd[variant_type][s]['frequencies']

            for samples_s_i in samples_s:

                if samples_s_i not in frequency_sample_dict:
                    frequency_sample_dict[samples_s_i] = {}
                    frequency_sample_dict[samples_s_i]['sites'] = []
                    frequency_sample_dict[samples_s_i]['frequencies'] = []

                # add the SITE
                frequency_sample_dict[samples_s_i]['sites'].append(s)
                ###
                if samples_s_i in samples_non_zero_s:
                    frequency_sample_dict[samples_s_i]['frequencies'].append(frequencies_s[samples_non_zero_s.index(samples_s_i)])
                else:
                    frequency_sample_dict[samples_s_i]['frequencies'].append(0)

        samples_all = list(frequency_sample_dict.keys())
        # turn freqs into numpy array to speed things up
        frequency_sample_dict_keys = list(frequency_sample_dict.keys())
        for s in frequency_sample_dict_keys:
            frequency_sample_dict[s]['sites'] = numpy.asarray(frequency_sample_dict[s]['sites'])
            frequency_sample_dict[s]['frequencies'] = numpy.asarray(frequency_sample_dict[s]['frequencies'])

        samples_pairs = itertools.combinations(samples_all, 2)
        error_range = numpy.logspace(numpy.log10(min(error_to_keep)), numpy.log10(max(error_to_keep)), num=20, base=10.0)

        genetic_diversity_dict = {}
        error_range_to_keep = []
        mean_pi_all = []
        mean_pi_ci_lower_all = []
        mean_pi_ci_upper_all = []
        percentile_all = []
        mean_n_all = []

        error_range = [0.01]

        for error_i in error_range:

            site_mapgd_to_keep_i = site_mapgd_to_keep[error_to_keep>=error_i]

            # need at least 100 sites across all samples
            if len(site_mapgd_to_keep_i) < min_n_sites:
                continue

            for samples_pair in samples_pairs:

                sample_1, sample_2 = samples_pair

                sample_1_sites = frequency_sample_dict[sample_1]['sites']
                sample_2_sites = frequency_sample_dict[sample_2]['sites']

                idx_sample_ = numpy.in1d(sample_i_sites, site_mapgd_to_keep_i)



                #    sample_i_frequencies = frequency_sample_dict[sample_i]['frequencies']



            #pi_all = []
            #n_all = []
            #for sample_i in samples_all:

            #    sample_i_sites = frequency_sample_dict[sample_i]['sites']
            #    sample_i_frequencies = frequency_sample_dict[sample_i]['frequencies']

            #    idx_sample_i = numpy.in1d(sample_i_sites, site_mapgd_to_keep_i)
            #    #idx_sample_i = numpy.asarray([sample_i_sites.index(s_) for s_ in site_mapgd_to_keep if s_ in sample_i_sites])
            #    sample_i_frequencies = sample_i_frequencies[idx_sample_i]

            #    #if len(sample_i_frequencies) == 0:
            #    #    pi = 0
            #    #else:
            #    pi = 2*sum(sample_i_frequencies*(1-sample_i_frequencies))/len(sample_i_frequencies)
            #    pi_all.append(pi)
            #    n_all.append(len(sample_i_frequencies))


            #mean_pi_null_all = []
            #for i in range(iter):
            #    # permute errors
            #    error_to_keep_null = numpy.random.permutation(error_to_keep)
            #    site_mapgd_to_keep_i_null = site_mapgd_to_keep[error_to_keep_null>=error_i]

            #    pi_null_all = []
            #    for sample_i in samples_all:

            #        sample_i_sites = frequency_sample_dict[sample_i]['sites']
            #        sample_i_frequencies = frequency_sample_dict[sample_i]['frequencies']
            #        idx_sample_i = numpy.in1d(sample_i_sites, site_mapgd_to_keep_i_null)
            #        sample_i_frequencies = sample_i_frequencies[idx_sample_i]
            #        pi_null = 2*sum(sample_i_frequencies*(1-sample_i_frequencies))/len(sample_i_frequencies)
            #        pi_null_all.append(pi_null)

            #    mean_pi_null_all.append(numpy.mean(pi_null_all))


            #mean_pi = numpy.mean(pi_all)
            #mean_n = numpy.mean(n_all)

            #mean_pi_null_all = numpy.asarray(mean_pi_null_all)
            #mean_pi_null_all = numpy.sort(mean_pi_null_all)

            #p_value = sum(slope_null_all > slope)/iter
            #lower_ci = mean_pi_null_all[int(iter*0.025)]
            #upper_ci = mean_pi_null_all[int(iter*0.975)]

            #percentile = sum(mean_pi_null_all < mean_pi)/iter

            #error_range_to_keep.append(error_i)
            #mean_pi_all.append(mean_pi)
            #mean_pi_ci_lower_all.append(lower_ci)
            #mean_pi_ci_upper_all.append(upper_ci)
            #percentile_all.append(percentile)
            #mean_n_all.append(mean_n)

            #print(error_i, mean_pi, percentile)
            #print(mean_pi_null_all)

        #error_vs_diversity_dict[variant_type]['error_range'] = error_range_to_keep
        #error_vs_diversity_dict[variant_type]['mean_pi'] = mean_pi_all
        #error_vs_diversity_dict[variant_type]['mean_pi_ci_lower'] = mean_pi_ci_lower_all
        #error_vs_diversity_dict[variant_type]['mean_pi_ci_upper'] = mean_pi_ci_upper_all
        #error_vs_diversity_dict[variant_type]['percentile'] = percentile_all
        #error_vs_diversity_dict[variant_type]['mean_n'] = mean_n_all

        #print(error_range_to_keep)
        #print(mean_pi_all)
        #print(percentile_all)



    error_vs_diversity_dict_path = error_vs_diversity_dict_template % species_name
    sys.stderr.write("Saving dict...\n")
    with open(error_vs_diversity_dict_path, 'wb') as handle:
        pickle.dump(error_vs_diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


        # Nei M 1986 Definition and esstimation of fixation indices. Evolution 40: 643-645
        # then, go through the error range and remove sites



species_list = ['Bacteroides_vulgatus_57955']

for species_name in species_list:

    make_error_vs_diversity_dict(species_name)
