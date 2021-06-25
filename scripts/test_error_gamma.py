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
from math import log10,ceil,fabs
from numpy.random import randint, choice, multinomial

import midas_db_utils

import scipy.stats as stats


import parse_HMP_data

allowed_variant_types = set(['1D','4D'])

good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = [good_species_list[3]]

import matplotlib.pyplot as plt

def generate_date():


    out_dict = {}


    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"


    for species_name in good_species_list:

        out_dict[species_name] = {}

        #intermediate_filename_template = config.data_directory+"pi_annotated_snps/%s.dat"

        #snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)
        snp_file_path = "%ssnps/%s/snps_summary.txt" % (config.data_directory, species_name)

        if os.path.isfile(snp_file_path) == False:
            continue

        samples = []

        for line in open(snp_file_path):

            line = line.strip().split('\t')
            sample = line[0]
            samples.append(sample)

        # Open post-processed MIDAS output
        #snp_file =  bz2.BZ2File(snp_file_path, "r")


        #line = snp_file.readline() # header
        #print(line)
        #items = line.split()[1:]
        #samples = numpy.array([item.strip() for item in items])

        samples = numpy.array(samples)

        #snp_file.close()

        #out_dict[species_name]['samples'] = samples


        richness = []
        diversity = []
        evenness = []
        #dominance = []

        # get strain data
        for sample in samples:

            intermediate_strain_filename = intermediate_strain_filename_template % sample

            if os.path.isfile(intermediate_strain_filename) == False:
                continue

            with open(intermediate_strain_filename, 'rb') as handle:
                b = pickle.load(handle)

            #print(species_name, b.keys())

            #print(b[species_name])

            if species_name in b:

                abundances = b[species_name]
                abundances = numpy.asarray(abundances)

                richness.append(len(abundances))

                if len(abundances) > 1:

                    H = -1*sum(abundances*numpy.log(abundances))
                    J = H/numpy.log(len(abundances))

                    diversity.append(H)
                    evenness.append(J)

        richness = numpy.asarray(richness)

        out_dict[species_name]['mean_richness'] = numpy.mean(richness)
        out_dict[species_name]['mean_diversity'] = numpy.mean(diversity)
        out_dict[species_name]['mean_evenness'] = numpy.mean(evenness)


        out_dict[species_name]['fraction_saples_strain'] = sum(richness>1)/len(richness)


        print(sum(richness>1)/len(richness), numpy.mean(richness),  numpy.mean(evenness))


        intermediate_filename_template = config.data_directory+"predicted_observed_prevalence/%s.dat"

        intermediate_filename = intermediate_filename_template % species_name

        if os.path.isfile(intermediate_filename) == False:
            continue


        with open(intermediate_filename, 'rb') as handle:
            snp_dict = pickle.load(handle)


        for allowed_variant in allowed_variant_types:

            out_dict[species_name][allowed_variant] = {}

            predicted_prevalence = snp_dict[allowed_variant]['predicted_prevalence']
            predicted_prevalence = numpy.asarray(predicted_prevalence)

            observed_prevalence = snp_dict[allowed_variant]['observed_prevalence']
            observed_prevalence = numpy.asarray(observed_prevalence)

            mean_abs_error = numpy.mean(numpy.absolute( observed_prevalence -  predicted_prevalence))

            mean_error = numpy.mean(observed_prevalence -  predicted_prevalence)

            out_dict[species_name][allowed_variant]['mean_abs_error'] = mean_abs_error
            out_dict[species_name][allowed_variant]['mean_error'] = mean_error


    output_filename = config.data_directory+"error_gamma.dat"

    with open(output_filename, 'wb') as handle:
        pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




output_filename = config.data_directory+"error_gamma.dat"

with open(output_filename, 'rb') as handle:
    error_dict = pickle.load(handle)


x = []
y_4D = []
y_1D = []
for species_name in good_species_list:

    if '4D' not in error_dict[species_name]:
        continue

    fraction_saples_strain = error_dict[species_name]['fraction_saples_strain']

    if (fraction_saples_strain==0) or (fraction_saples_strain == 1):
        continue


    mean_abs_error_4D = error_dict[species_name]['4D']['mean_abs_error']
    mean_abs_error_1D = error_dict[species_name]['1D']['mean_abs_error']

    print(species_name,fraction_saples_strain, mean_abs_error_4D)


    x.append(fraction_saples_strain)
    y_4D.append(mean_abs_error_4D)
    y_1D.append(mean_abs_error_1D)




fig, ax = plt.subplots(figsize=(4,4))

ax.scatter(x, y_4D, c='b', alpha=0.7, label='Species, 4D')
#ax.scatter(x, y_1D, c='r', alpha=0.7, label='Species, 1D')


slope_4D, intercept_4D, r_value_4D, p_value_4D, std_err_4D = stats.linregress(numpy.log10(x), numpy.log10(y_4D))
slope_1D, intercept_1D, r_value_1D, p_value_1D, std_err_1D = stats.linregress(numpy.log10(x), numpy.log10(y_1D))


print(p_value_4D, p_value_1D)

x_range =  numpy.linspace(min(numpy.log10(x)), max(numpy.log10(x)), 10000)

y_fit_range_4D = (slope_4D*x_range + intercept_4D)
y_fit_range_1D = (slope_1D*x_range + intercept_1D)

ax.plot(10**x_range, 10**y_fit_range_4D, c='b', lw=2.5, linestyle='--', zorder=2)
#ax.plot(10**x_range, 10**y_fit_range_1D, c='r', lw=2.5, linestyle='--', zorder=2)


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel("Fraction hosts with > 1 strain", fontsize = 12)
ax.set_ylabel("Mean absolute error", fontsize = 12)


ax.legend(loc="upper left", fontsize=8)


fig.tight_layout()
fig.savefig("%s%s.png" % (config.analysis_directory, 'strain_error_gamma'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



#print(error_dict)
