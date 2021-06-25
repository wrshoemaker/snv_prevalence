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


random.seed(123456789)


#species_name='Bacteroides_vulgatus_57955'

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = ['Bacteroides_vulgatus_57955']


def jaccard(x,y):
  x = numpy.asarray(x, numpy.bool) # Not necessary, if you keep your data
  y = numpy.asarray(y, numpy.bool) # in a boolean array already!
  return numpy.double(numpy.bitwise_and(x, y).sum()) / numpy.double(numpy.bitwise_or(x, y).sum())


for species_name in good_species_list:

    sys.stderr.write("%s\n" % species_name)

    ref_freq_file_path = "%ssnps/%s/snps_ref_freq.txt.bz2" % (data_directory, species_name)
    alt_allele_file_path = "%ssnps/%s/snps_alt_allele.txt.bz2" % (data_directory, species_name)
    info_file_path = "%ssnps/%s/snps_info.txt.bz2" % (data_directory, species_name)

    if os.path.isfile(ref_freq_file_path) == False:
        continue

    if os.path.isfile(alt_allele_file_path) == False:
        continue

    if os.path.isfile(info_file_path) == False:
        continue

    ref_freq_file = bz2.BZ2File(ref_freq_file_path,"r")
    alt_allele_file = bz2.BZ2File(alt_allele_file_path,"r")
    info_file = bz2.BZ2File(info_file_path,"r")

    #ref_freq_line = ref_freq_file.readline().decode('utf-8')
    #alt_line = alt_allele_file.readline().decode('utf-8')
    #info_line = info_file.readline().decode('utf-8')

    ref_freq_line = ref_freq_file.readline()
    alt_line = alt_allele_file.readline()

    alt_line_items = alt_line.split()
    samples = numpy.array(alt_line_items[1:])
    samples_set = set(samples)

    # filter samples so no twins, etc


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



    #sample_pair_freqs_dict = {}

    #for allowed_variant in allowed_variant_types:
    #    sample_pair_freqs_dict[allowed_variant] = {}

    #for sample_pair in sample_pairs:
    #    key_name = '%s_%s' % sample_pair
    #    for allowed_variant in allowed_variant_types:
    #        sample_pair_freqs_dict[allowed_variant][key_name] = {}
    #        sample_pair_freqs_dict[allowed_variant][key_name]['sample_1'] = []
    #        sample_pair_freqs_dict[allowed_variant][key_name]['sample_2'] = []

    sys.stderr.write("Creating 2D SFS dictionary...\n")
    count = 0

    site_occupancy_dict = {}
    for allowed_variant_type in allowed_variant_types:
        site_occupancy_dict[allowed_variant_type] = {}

    for ref_freq_line in ref_freq_file:

        if (count % 100000 == 0) and (count > 0):
            sys.stderr.write("%d sites: complete\n" % count)

        #if (count % 100 == 0) and (count > 0):
        #    sys.stderr.write("%d sites: complete\n" % count)

        #if count >= 10000:
        #    break

        count += 1

        if count > 5000:
            continue



        #ref_freq_line = ref_freq_file.readline().decode('utf-8')
        #ref_freq_line = ref_freq_line.decode('utf-8')
        ref_freq_position = ref_freq_line.split()[0]
        ref_freq_position_aa = site_aa_dict[ref_freq_position]

        if ref_freq_position_aa not in allowed_variant_types:
            continue

        ref_freqs = numpy.array([float(item) for item in ref_freq_line.split()[1:]])

        # go through and find the samples that are polymorphic to reduce runtime
        # only loop through ~10**2 pairs instead of ~1000**2 pairs
        line_sample_freq_dict = {sample:ref_freqs[sample_idx] for sample_idx, sample in enumerate(samples) if (ref_freqs[sample_idx]>0) and (ref_freqs[sample_idx]<1)}
        #line_sample_freq_dict = {sample:ref_freqs[sample_idx] for sample_idx, sample in enumerate(samples) }
        #print(len(line_sample_freq_dict), len(samples))
        if len(line_sample_freq_dict) == 0:
            continue

        for sample_ in line_sample_freq_dict.keys():
            if sample_ not in samples:
                line_sample_freq_dict.pop(sample_, None)

            line_sample_freq_dict[sample_] = 1

        #for sample_ in samples:
        #    if

        samples_diff = list(samples_set - set(line_sample_freq_dict.keys()))
        for sample_ in samples_diff:
            line_sample_freq_dict[sample_] = 0

        site_occupancy_dict[ref_freq_position_aa][ref_freq_position] = line_sample_freq_dict

    sys.stderr.write("next line...\n")


    distances_dict = {}
    for allowed_variant_type in allowed_variant_types:
        distances_dict[allowed_variant_type] = {}
        distances_dict[allowed_variant_type]['site_occupancies'] = {}
        distances_dict[allowed_variant_type]['site_pairs'] = []
        distances_dict[allowed_variant_type]['jaccard_distances'] = []
        distances_dict[allowed_variant_type]['pearsons_rho'] = []

    sys.stderr.write("dictioanry init...\n")

    for allowed_variant_type in allowed_variant_types:
        #site_pairs = list(itertools.combinations(list(site_occupancy_dict[allowed_variant_type].keys()), 2))
        sites = list(site_occupancy_dict[allowed_variant_type].keys())
        for site_i_idx, site_i in enumerate(sites):
            #print(allowed_variant_type, site_i_idx)
            #sys.stderr.write("%s %s...\n" % (str(allowed_variant_type), str(site_i_idx)))
            # get occupanceies
            site_i_occupancies = [ site_occupancy_dict[allowed_variant_type][site_i][sample] for sample in samples ]
            site_i_occupancies = numpy.asarray(site_i_occupancies)
            site_i_occupancy = sum(site_i_occupancies) / len(site_i_occupancies)

            if site_i_occupancy > float(0):
                distances_dict[allowed_variant_type]['site_occupancies'][site_i] = site_i_occupancy

                #print(site_i, site_i_occupancy)
                for site_j_idx, site_j in enumerate(sites):

                    if site_j_idx <= site_i_idx:
                        continue

                    site_j_occupancies = [ site_occupancy_dict[allowed_variant_type][site_j][sample] for sample in samples ]
                    site_j_occupancies = numpy.asarray(site_j_occupancies)

                    distance_ij = jaccard(site_i_occupancies, site_j_occupancies)
                    rho_ij = numpy.corrcoef(site_i_occupancies, site_j_occupancies)[0,1]

                    distances_dict[allowed_variant_type]['site_pairs'].append('%s_%s' % (site_i, site_j))
                    distances_dict[allowed_variant_type]['jaccard_distances'].append(distance_ij)
                    distances_dict[allowed_variant_type]['pearsons_rho'].append(rho_ij)



    ref_freq_file.close()
    alt_allele_file.close()
    info_file.close()

    sys.stderr.write("Done!\n")

    sys.stderr.write("Writing distance file...\n")


    filename_ = "%sjaccard_distance_snps/%s.pickle" % (data_directory, species_name)

    with open(filename_, 'wb') as handle:
        pickle.dump(distances_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
