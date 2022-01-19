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


count = 100

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

#color_dict = {'4D': 'b', '1D': 'r'}

intermediate_filename_template = config.data_directory+"coprevalence_f0/%s.dat"


good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = ['Bacteroides_vulgatus_57955']

f0_range = numpy.logspace(numpy.log10(0.01), numpy.log10(0.2), num=20, endpoint=True, base=10.0)


for species_name in good_species_list:

    # Holds panel wide prevalence for each species
    os.system('mkdir -p %ssnp_prevalences' % config.data_directory)

    snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

    if os.path.isfile(snp_file_path) == False:
        continue

    # Open post-processed MIDAS output
    snp_file =  bz2.BZ2File(snp_file_path, "r")


    line = snp_file.readline() # header
    items = line.split()[1:]
    samples = numpy.array([item.strip() for item in items])

    count_variant_type_dict = {}
    site_occupancy_dict = {}
    for allowed_variant_type in allowed_variant_types:
        site_occupancy_dict[allowed_variant_type] = {}
        count_variant_type_dict[allowed_variant_type] = 0


    sys.stderr.write("Calculating for %s\n" % species_name)


    sys.stderr.write("Calculating SNP prevalences...\n")
    num_sites_processed = 0
    for line in snp_file:

        num_sites_processed+=1
        #print num_sites_processed
        if num_sites_processed%50000==0:
            sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))
            #if debug:
            #    break

        items = line.split()
        #print(items)
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = long(info_items[1])
        gene_name = info_items[2]
        location_aa = info_items[3]

        if location_aa not in allowed_variant_types:
            continue


        if (count_variant_type_dict['1D'] >= count) and (count_variant_type_dict['4D'] >= count):
            break

        if (count_variant_type_dict[location_aa] >= count):
            continue



        #while (count_variant_type_dict['1D'] < count) and (count_variant_type_dict['4D'] < count):
        #print(location_aa, count_variant_type_dict[location_aa])

        #
        # Load alt and depth counts
        alts = []
        depths = []
        for item in items[1:]:
            subitems = item.split(",")
            alts.append(long(subitems[0]))
            depths.append(long(subitems[1]))
        alts = numpy.array(alts)
        depths = numpy.array(depths)
        refs = depths-alts

        # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
        population_prevalence = ((alts>=refs)*(depths>0)).sum()
        population_freq = population_prevalence*1.0/(depths>0).sum()

        if population_freq>0.5:
            # alternate allele is in the majority
            # re-polarize for now
            alts,refs = refs,alts

        # Next calculate fraction of cohort where population minor allele is present at >=10% within-host frequency
        alt_threshold = numpy.ceil(depths*0.1)+0.5 #at least one read above 10%.

        snp_prevalence = ((alts>=alt_threshold)*(depths>0)).sum()
        snp_freq = snp_prevalence*1.0/(depths>0).sum()

        if (population_prevalence==0) and (snp_prevalence==0):
            continue

        if (population_prevalence==1) and (snp_prevalence==1):
            continue

        if (snp_prevalence==0) or (snp_prevalence==1):
            continue

        site_occupancy_dict[location_aa][location] = (alts>=alt_threshold)*(depths>0)

        count_variant_type_dict[location_aa] += 1




    #continue
    #fig, ax = plt.subplots(figsize=(4,4))

    pairwise_dict = {}

    for allowed_variant_type in allowed_variant_types:


        #ratios = []

        pairwise_dict[allowed_variant_type] = {}

        for f0 in f0_range:

            #print(allowed_variant_type, f0)

            #sys.stderr.write("Calculating for %s\n" % species_name)


            #site_pairs = list(itertools.combinations(list(site_occupancy_dict[allowed_variant_type].keys()), 2))
            locations = list(site_occupancy_dict[allowed_variant_type].keys())

            #correlation_list = []
            #minimum_prevalence_list = []

            sigma_k_numerator_list = []
            sigma_k_denominator_list = []

            for location_i_idx, location_i in enumerate(locations):
                #sys.stderr.write("%s %s...\n" % (str(allowed_variant_type), str(site_i_idx)))
                # get occupanceies
                #site_i_occupancies = [ site_occupancy_dict[allowed_variant_type][site_i][sample] for sample in samples ]
                #site_i_occupancies = numpy.asarray(site_i_occupancies)
                #site_i_occupancy = sum(site_i_occupancies) / len(site_i_occupancies)

                prevalence_array_i = site_occupancy_dict[allowed_variant_type][location_i]
                prevalence_i = sum(prevalence_array_i) / len(prevalence_array_i)

                #print(site_i, site_i_occupancy)
                for location_j_idx, location_j in enumerate(locations):

                    if location_j_idx <= location_i_idx:
                        continue

                    prevalence_array_j = site_occupancy_dict[allowed_variant_type][location_j]
                    prevalence_j = sum(prevalence_array_j) / len(prevalence_array_j)


                    #joint_prevalence_array = numpy.logical_and(prevalence_array_i, prevalence_array_j)
                    #joint_prevalence = sum(joint_prevalence_array) / len(joint_prevalence_array)
                    #numerator = (joint_prevalence - prevalence_i*prevalence_j)**2
                    #denominator = prevalence_i*(1-prevalence_i)*prevalence_j*(1-prevalence_j)
                    #correlation = numerator/denominator
                    #minimum_prevalence = min([prevalence_i, prevalence_j])
                    #if denominator <= 0:
                    #    print(minimum_prevalence, denominator)
                    #correlation_list.append(correlation)
                    #minimum_prevalence_list.append(minimum_prevalence)

                    sigma_k_numerator, sigma_k_denominator = diversity_utils.calculate_unbiased_sigma_k_vector(prevalence_array_i, prevalence_array_j, f0, k=2)

                    if (math.isnan(sigma_k_numerator) == True) or (math.isnan(sigma_k_denominator) == True):
                        continue

                    #print(sigma_k_numerator, sigma_k_denominator)

                    sigma_k_numerator_list.append(sigma_k_numerator)
                    sigma_k_denominator_list.append(sigma_k_denominator)


            #sigma_k_numerator_list = numpy.asarray(sigma_k_numerator_list)
            #sigma_k_denominator_list = numpy.asarray(sigma_k_denominator_list)

            # output numerator and denominator as dict

            pairwise_dict[allowed_variant_type][f0] = {}
            pairwise_dict[allowed_variant_type][f0]['sigma_k_numerator'] = sigma_k_numerator_list
            pairwise_dict[allowed_variant_type][f0]['sigma_k_denominator'] = sigma_k_denominator_list



            #ratios.append(sum(sigma_k_numerator_list) / sum(sigma_k_denominator_list))
        #ratios = numpy.asarray(ratios)
        #ax.scatter(f0_range, ratios, c=color_dict[allowed_variant_type])

    intermediate_filename = intermediate_filename_template % species_name

    with open(intermediate_filename, 'wb') as handle:
        pickle.dump(pairwise_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



    #ax.set_xscale('log', basex=10)
    #ax.set_yscale('log', basey=10)

    #ax.set_xlabel("Weighting frequency, " + r'$f_{0}$', fontsize = 12)
    #ax.set_ylabel("SNP prevalence disequilibria, " + r'$\sigma^{2}_{d}$', fontsize = 12)

    #ax.set_title("SNP prevalence")

    #ax.set_xlim([min(f0_range)*0.8, max(f0_range)*1.2])
    #ax.set_ylim([min(ratios)*0.8, max(ratios)*1.2])

    #fig.tight_layout()
    #fig.savefig("%s%s.png" % (config.analysis_directory, 'prevalence_correlation_vs_f0'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    #plt.close()
