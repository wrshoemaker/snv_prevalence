from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import parse_midas_data

import diversity_utils

import matplotlib.pyplot as plt


data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])


good_species_list = parse_midas_data.parse_good_species_list()

color_dict = {'4D': 'b', '1D': 'r'}


intermediate_filename_template = config.data_directory+"coprevalence_f0/%s.dat"


for species_name in good_species_list:

    intermediate_filename = intermediate_filename_template % species_name

    with open(intermediate_filename, 'rb') as handle:
        coprevalence_dict = pickle.load(handle)



    is_nan = False

    dicttt = {}

    for allowed_variant_type in allowed_variant_types:

        f0_list = list(coprevalence_dict[allowed_variant_type].keys())
        f0_list.sort()
        f0_list = numpy.asarray(f0_list)

        f0_numerator = [sum(coprevalence_dict[allowed_variant_type][f0]['sigma_k_numerator']) for f0 in f0_list]
        f0_denominator = [sum(coprevalence_dict[allowed_variant_type][f0]['sigma_k_denominator']) for f0 in f0_list]

        f0_numerator = numpy.asarray(f0_numerator)
        f0_denominator = numpy.asarray(f0_denominator)

        ratios = f0_numerator/f0_denominator

        if (sum(numpy.isnan(ratios)) > 0) :
            is_nan = True




        dicttt[allowed_variant_type] = {}
        dicttt[allowed_variant_type]['f0'] = f0_list
        dicttt[allowed_variant_type]['ratios'] = ratios


    if is_nan == True:
        continue

    fig, ax = plt.subplots(figsize=(4,4))

    f0_max = []
    ratios_max = []

    for allowed_variant_type in allowed_variant_types:
        ratios = dicttt[allowed_variant_type]['ratios']
        f0_list = dicttt[allowed_variant_type]['f0']
        ax.scatter(f0_list, ratios, c=color_dict[allowed_variant_type])

        f0_max.append(max(f0_list))
        f0_max.append(min(f0_list))

        ratios_max.append(max(ratios))
        ratios_max.append(min(ratios))

    print(f0_max)
    print(ratios_max)

    ax.set_xlim([min(f0_max)*0.8, max(f0_max)*1.2])
    ax.set_ylim([min(ratios_max)*0.8, max(ratios_max)*1.2])

    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

    ax.set_xlabel("Weighting frequency, " + r'$f_{0}$', fontsize = 12)
    ax.set_ylabel("SNP prevalence disequilibria, " + r'$\sigma^{2}_{d}$', fontsize = 12)

    ax.set_title(species_name)



    fig.tight_layout()
    fig.savefig("%scoprevalence_f0/%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
