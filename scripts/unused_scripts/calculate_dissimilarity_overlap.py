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
min_muts=10
#samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D']))

#sample_pair_example = ('700013597', '700015702')

#good_species_list = parse_midas_data.parse_good_species_list()

good_species_list = ['Bacteroides_vulgatus_57955']

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

    sample_pairs = list(itertools.combinations(samples, 2))

    #sample_pairs = [random.choice(sample_pairs)]

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



    sample_pair_freqs_dict = {}

    for allowed_variant in allowed_variant_types:
        sample_pair_freqs_dict[allowed_variant] = {}

    for sample_pair in sample_pairs:
        key_name = '%s_%s' % sample_pair
        for allowed_variant in allowed_variant_types:
            sample_pair_freqs_dict[allowed_variant][key_name] = {}
            sample_pair_freqs_dict[allowed_variant][key_name]['sample_1'] = []
            sample_pair_freqs_dict[allowed_variant][key_name]['sample_2'] = []

    sys.stderr.write("Creating 2D SFS dictionary...\n")
    count = 0


    for ref_freq_line in ref_freq_file:

        if (count % 100000 == 0) and (count > 0):
            sys.stderr.write("%d sites: complete\n" % count)

        #if count > 10000:
        #    break

        count += 1

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

        line_sample_pairs = list(itertools.combinations(list(line_sample_freq_dict.keys()), 2))


        for sample_pair in line_sample_pairs:

            key_name = '%s_%s' % sample_pair

            if key_name not in sample_pair_freqs_dict[ref_freq_position_aa]:
                continue

            freq_sample_1 = line_sample_freq_dict[sample_pair[0]]
            freq_sample_2 = line_sample_freq_dict[sample_pair[1]]

            if freq_sample_1 < 0.5:
                freq_sample_1 = 1 - freq_sample_1

            if freq_sample_2 < 0.5:
                freq_sample_2 = 1 - freq_sample_2


            sample_pair_freqs_dict[ref_freq_position_aa][key_name]['sample_1'].append(freq_sample_1)
            sample_pair_freqs_dict[ref_freq_position_aa][key_name]['sample_2'].append(freq_sample_2)



    #print(sample_pair_freqs_dict)

    ref_freq_file.close()
    alt_allele_file.close()
    info_file.close()

    sys.stderr.write("Done!\n")

    sys.stderr.write("Writing deviation file...\n")


    filename_sfs_2D = "%sdissimilarity_overlap/%s.txt.gz" % (data_directory, species_name)

    record_strs = [", ".join(['variant_type', 'sample_1', 'sample_2', 'overlap', 'D_rJSD'])]


    for variant_type, variant_dict in sample_pair_freqs_dict.items():

        for sample_pair in variant_dict.keys():

            sample_1 = sample_pair.split('_')[0]
            sample_2 = sample_pair.split('_')[1]

            freqs_1 = variant_dict[sample_pair]['sample_1']
            freqs_2 = variant_dict[sample_pair]['sample_2']

            freqs_1 = numpy.asarray(freqs_1)
            freqs_2 = numpy.asarray(freqs_2)

            if (len(freqs_1) == 0) or (len(freqs_1) == 0):
                continue

            freqs_1_shared = freqs_1[(freqs_1<1) & (freqs_2<1)]
            freqs_2_shared = freqs_2[(freqs_1<1) & (freqs_2<1)]

            if len(freqs_2_shared) < min_muts:
                continue

            # divide by number sites so overlap is between zero and 1
            overlap = sum(freqs_1_shared + freqs_2_shared)/(2 * len(freqs_1_shared))

            freqs_1_shared_renorm = freqs_1_shared / sum(freqs_1_shared)
            freqs_2_shared_renorm = freqs_2_shared / sum(freqs_2_shared)

            m = (freqs_1_shared_renorm + freqs_2_shared_renorm)/2

            D_kl_1 = sum(freqs_1_shared_renorm*numpy.log(freqs_1_shared_renorm/m))
            D_kl_2 = sum(freqs_2_shared_renorm*numpy.log(freqs_2_shared_renorm/m))

            D_rJSD = ((D_kl_1 + D_kl_2)/2) ** (1/2)

            record_sample_pair_str = ", ".join([variant_type, sample_1, sample_2, str(overlap), str(D_rJSD)])

            record_strs.append(record_sample_pair_str)


    file_sfs_2d = gzip.open(filename_sfs_2D,"w")

    record_str = "\n".join(record_strs)

    file_sfs_2d.write(record_str.encode())
    file_sfs_2d.close()

    sys.stderr.write("Done!\n")
