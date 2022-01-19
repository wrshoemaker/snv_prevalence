from __future__ import division
import os, sys
import bz2
import random
import itertools
import config
import parse_midas_data
import numpy
import pickle

import gzip

import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import gamma

random.seed(123456789)


min_coverage = 50
count = 1000

min_frequency = 0.025

#species_name='Bacteroides_vulgatus_57955'

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

#good_species_list = parse_midas_data.parse_good_species_list()

good_species_list = ['Bacteroides_vulgatus_57955']


for species_name in good_species_list:

    sys.stderr.write("%s\n" % species_name)

    ref_freq_file_path = "%ssnps/%s/snps_ref_freq.txt.bz2" % (data_directory, species_name)
    alt_allele_file_path = "%ssnps/%s/snps_alt_allele.txt.bz2" % (data_directory, species_name)
    info_file_path = "%ssnps/%s/snps_info.txt.bz2" % (data_directory, species_name)
    depth_file_path = "%ssnps/%s/snps_depth.txt.bz2" % (data_directory, species_name)


    if os.path.isfile(ref_freq_file_path) == False:
        continue

    if os.path.isfile(alt_allele_file_path) == False:
        continue

    if os.path.isfile(info_file_path) == False:
        continue

    if os.path.isfile(depth_file_path) == False:
        continue


    ref_freq_file = bz2.BZ2File(ref_freq_file_path,"r")
    alt_allele_file = bz2.BZ2File(alt_allele_file_path,"r")
    info_file = bz2.BZ2File(info_file_path,"r")
    depth_file = bz2.BZ2File(depth_file_path,"r")

    ref_freq_line = ref_freq_file.readline()
    alt_line = alt_allele_file.readline()
    info_line = info_file.readline()
    depth_line = depth_file.readline()

    alt_line_items = alt_line.split()
    samples = numpy.array(alt_line_items[1:])
    samples_set = set(samples)

    samples_idx_dict = {sample:sample_idx for sample_idx, sample in enumerate(samples) }

    site_aa_dict = {}
    sys.stderr.write("Creating amino acid redundancy map...\n")

    for info_file_line in info_file:
        #info_file_line = info_file.readline().decode('utf-8')
        #info_file_line = info_file_line.decode('utf-8')
        info_file_items = info_file_line.split()

        if len(info_file_items) == 0:
            continue

        if info_file_items[-1] == 'NC':
            site_aa_dict[info_file_items[0]] = info_file_items[-1]

        elif ('SYN' in info_file_items[-1]) or ('NS' in info_file_items[-1]):
            site_aa_dict[info_file_items[0]] = info_file_items[5]

    sys.stderr.write("Done!\n")


    count_variant_dict = {}
    site_occupancy_dict = {}
    for allowed_variant_type in allowed_variant_types:
        site_occupancy_dict[allowed_variant_type] = {}

        count_variant_dict[allowed_variant_type] = 0

    # create dictionary with coverage
    position_coverage_dict = {}
    for depth_line_idx, depth_line in enumerate(depth_file):

        if depth_line_idx > 1000000:
            break

        depth_position = depth_line.split()[0]
        depth_position_aa = site_aa_dict[depth_position]

        if depth_position_aa not in allowed_variant_types:
            continue

        if count_variant_dict[depth_position_aa] > count:
            continue

        count_variant_dict[depth_position_aa] += 1

        depth_values = numpy.array([float(item) for item in depth_line.split()[1:]])

        #depth_values_keep = [depth_values>min_coverage]

        position_coverage_dict[depth_position] = depth_values


    freqs = []

    mean_freqs_no_zeros_list = []
    var_freqs_no_zeros_list = []
    prevalence_list = []

    for ref_freq_line_idx, ref_freq_line in enumerate(ref_freq_file):

        if ref_freq_line_idx > 1000000:
            break

        ref_freq_position = ref_freq_line.split()[0]
        ref_freq_position_aa = site_aa_dict[ref_freq_position]

        if ref_freq_position_aa not in allowed_variant_types:
            continue

        if ref_freq_position not in position_coverage_dict:
            continue

        ref_freqs = numpy.array([float(item) for item in ref_freq_line.split()[1:]])

        ref_freqs_depths = position_coverage_dict[ref_freq_position]

        ref_count_alt = ref_freqs*ref_freqs_depths
        ref_count_alt = ref_count_alt.astype(numpy.int)
        ref_freqs_depths = ref_freqs_depths.astype(numpy.int)

        ref_freqs_depths_filter = ref_freqs_depths[ref_freqs_depths >= min_coverage]
        ref_count_alt_filter = ref_count_alt[ref_freqs_depths >= min_coverage]
        ref_freqs_filter = ref_freqs[ref_freqs_depths >= min_coverage]


        # fold the frequencies
        ref_freqs_filter_greater = numpy.where(ref_freqs_filter >= 0.5)
        ref_freqs_filter[ref_freqs_filter_greater] = 1 - ref_freqs_filter[ref_freqs_filter_greater]

        ref_freqs_filter = ref_freqs_filter[ref_freqs_filter>0]
        ref_freqs_filter_list = ref_freqs_filter.tolist()

        if ref_freq_position_aa == '4D':

            freqs.extend(ref_freqs_filter_list)

            freqs_no_zeros_idx = ref_freqs_filter>0
            #freqs_no_zeros_idx = ref_freqs_filter>=0


            if sum(freqs_no_zeros_idx) > 3:

                freqs_no_zeros = ref_freqs_filter[freqs_no_zeros_idx]

                mean_freq_no_zero = numpy.mean(freqs_no_zeros)
                mean_freqs_no_zeros_list.append(mean_freq_no_zero)

                var_freqs_no_zeros_list.append(numpy.var( ref_freqs_filter[freqs_no_zeros_idx] ))


                #prevalence = sum(freqs_no_zeros_idx) / len(freqs_no_zeros_idx)


                #prevalence = sum(ref_freqs_filter>0) / len(ref_freqs_filter)

                prevalence = len(freqs_no_zeros) / len(ref_freqs_filter)
                print(prevalence)
                prevalence_list.append(prevalence)

                # min_frequency



        #ref_freqs_to_keep = ref_freqs[ref_freq_depths>min_coverage]


        #i = np.where(current_img_list == 82459)
        #current_img_list[i] -= 1

        #i = numpy.where(x >= 0.5)
        #x[i] = 1 - x[i]


        #run_sites = sum(count_variant_dict.values())


        #ref_freqs = numpy.array([float(item) for item in ref_freq_line.split()[1:]])
        #ref_freqs_minor = ref_freqs[(ref_freqs<0.5) & (ref_freqs>0.0)]
        #ref_freqs_minor_list = ref_freqs_minor.tolist()

        #ref_freqs_major = ref_freqs[(ref_freqs>0.5) & (ref_freqs<1.0)]
        #ref_freqs_major = 1-ref_freqs_major
        #ref_freqs_major_list = ref_freqs_major.tolist()

        #freqs.extend(ref_freqs_minor_list)
        #freqs.extend(ref_freqs_major_list)


    ref_freq_file.close()
    alt_allele_file.close()
    info_file.close()
    depth_file.close()

    sys.stderr.write("Done!\n")

    freqs = numpy.asarray(freqs)

    mean_freqs_no_zeros_list = numpy.asarray(mean_freqs_no_zeros_list)
    var_freqs_no_zeros_list = numpy.asarray(var_freqs_no_zeros_list)
    prevalence_list = numpy.asarray(prevalence_list)

    #freqs = freqs[(freqs>0.1) & (freqs<0.4)]

    freqs = freqs[freqs>0.07]

    freqs_log10 = numpy.log10(freqs)

    ag,bg,cg = gamma.fit(freqs_log10)

    x_range = numpy.linspace(min(freqs_log10) , max(freqs_log10) , 10000)





    # regressions


    mean_freqs_no_zeros_list_1 = mean_freqs_no_zeros_list[mean_freqs_no_zeros_list < min_frequency]
    var_freqs_no_zeros_list_1 = var_freqs_no_zeros_list[mean_freqs_no_zeros_list < min_frequency]

    mean_freqs_no_zeros_list_2 = mean_freqs_no_zeros_list[mean_freqs_no_zeros_list >= min_frequency]
    var_freqs_no_zeros_list_2 = var_freqs_no_zeros_list[mean_freqs_no_zeros_list >= min_frequency]

    prevalences_2 = prevalence_list[mean_freqs_no_zeros_list >= min_frequency]

    #print(prevalences_2)

    slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(numpy.log10(mean_freqs_no_zeros_list_1), numpy.log10(var_freqs_no_zeros_list_1))
    slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(numpy.log10(mean_freqs_no_zeros_list_2), numpy.log10(var_freqs_no_zeros_list_2))


    fig, ax = plt.subplots(figsize=(4,4))

    #ax.hist(freqs_log10, bins=30, density=True, alpha=0.7)
    #ax.plot(x_range, gamma.pdf(x_range, ag, bg,cg), 'k', lw=2)


    ax.scatter(mean_freqs_no_zeros_list, var_freqs_no_zeros_list, alpha=0.2)


    x_log10_range_1 =  numpy.linspace(min(numpy.log10(mean_freqs_no_zeros_list_1)) , max(numpy.log10(mean_freqs_no_zeros_list_1)) , 10000)
    y_log10_fit_range_1 = 10 ** (slope_1 * x_log10_range_1 + intercept_1)
    ax.plot(10**x_log10_range_1, y_log10_fit_range_1, c='k', lw=2.5, linestyle='--', zorder=2, label="Slope=%s" % str(round(slope_1, 3)))

    x_log10_range_2 =  numpy.linspace(min(numpy.log10(mean_freqs_no_zeros_list_2)) , max(numpy.log10(mean_freqs_no_zeros_list_2)) , 10000)
    y_log10_fit_range_2 = 10 ** (slope_2 * x_log10_range_2 + intercept_2)
    ax.plot(10**x_log10_range_2, y_log10_fit_range_2, c='k', lw=2.5, linestyle=':', zorder=2, label="Slope=%s" % str(round(slope_2, 3)))

    #ax.axvline(x=0.02)
    ax.axvline(x=min_frequency)
    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

    ax.set_xlim([min(mean_freqs_no_zeros_list)*0.8, max(mean_freqs_no_zeros_list)*1.2])
    ax.set_ylim([min(var_freqs_no_zeros_list)*0.8, max(var_freqs_no_zeros_list)*1.2])

    ax.set_xlabel("Mean SNP frequency, no zeros", fontsize = 12)
    ax.set_ylabel("Variance of SNP frequency, no zeros", fontsize = 12)

    ax.set_title(species_name)

    ax.legend(loc="lower right", fontsize=8)


    fig.tight_layout()
    fig.savefig("%s%s.png" % (config.analysis_directory, 'mean_vs_variance'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





    fig, ax = plt.subplots(figsize=(4,4))


    ax.scatter(mean_freqs_no_zeros_list, prevalence_list)


    fig.tight_layout()
    fig.savefig("%s%s.png" % (config.analysis_directory, 'mean_vs_prevalence'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
