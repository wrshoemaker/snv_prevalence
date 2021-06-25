from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import collections
import os.path
import scipy.stats as stats

import diversity_utils
import core_gene_utils
import parse_midas_data
import parse_HMP_data
import sample_utils
import calculate_substitution_rates
import clade_utils

import matplotlib.pyplot as plt
from scipy.stats import gamma


numpy.random.seed(123456789)

#count = 1000

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])
clade_types = ['all','largest_clade']

#color_dict = {'4D': 'b', '1D': 'r'}

#intermediate_filename_template = config.data_directory+"coprevalence_f0/%s.dat"
low_divergence_threshold = config.between_low_divergence_threshold

D_min=20

min_n_samples = 20
min_sample_size = config.between_host_min_sample_size


good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = ['Bacteroides_vulgatus_57955']


def get_samples_no_strains(species_name, samples):

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    samples_no_strains = []

    for sample in samples:

        intermediate_strain_filename = intermediate_strain_filename_template % sample

        if os.path.isfile(intermediate_strain_filename) == False:
            continue

        with open(intermediate_strain_filename, 'rb') as handle:
            b = pickle.load(handle)

        if species_name in b:

            abundances = b[species_name]

            if len(abundances) == 1:

                samples_no_strains.append(str(sample))

    samples_no_strains = numpy.asarray(samples_no_strains)

    return samples_no_strains




def calculate_pi_from_snps():

    #if include_samples_with_strains == False:
    #    intermediate_filename_template = config.data_directory+"pi_annotated_snps_no_strains/%s.dat"

    #else:
    intermediate_filename_template = config.data_directory+"pi_annotated_snps/%s.dat"

    for species_name in good_species_list:


        sys.stderr.write("Loading whitelisted genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)

        snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        if os.path.isfile(snp_file_path) == False:
            continue

        # Open post-processed MIDAS output
        snp_file =  bz2.BZ2File(snp_file_path, "r")
        line = snp_file.readline() # header
        items = line.split()[1:]
        samples = numpy.array([item.strip() for item in items])


        pi_dict = {}

        for allowed_variant_type in allowed_variant_types:
            pi_dict[allowed_variant_type] = {}

            for sample in samples:
                pi_dict[allowed_variant_type][sample] = {}

                pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] = 0
                pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] = 0

                pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] = 0
                pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] = 0


        sys.stderr.write("Iterating though SNPs...\n")
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

            # only examine core genes
            if gene_name not in core_genes:
                continue

            if location_aa not in allowed_variant_types:
                continue

            alts = []
            depths = []
            for item in items[1:]:
                subitems = item.split(",")
                alts.append(long(subitems[0]))
                depths.append(long(subitems[1]))
            alts = numpy.array(alts)
            depths = numpy.array(depths)
            refs = depths-alts


            # only try and predict prevalence if there's non-zero coverage for >= 20 hos
            if sum(depths>D_min) < min_n_samples:
                continue

            # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
            population_prevalence = ((alts>=refs)*(depths>D_min)).sum()
            population_freq = population_prevalence*1.0/(depths>D_min).sum()


            if population_freq>0.5:
                # alternate allele is in the majority
                # re-polarize for now
                alts,refs = refs,alts

            # Next calculate fraction of cohort where population minor allele is present at >=10% within-host frequency
            alt_threshold = numpy.ceil(depths*0.1)+0.5 #at least one read above 10%.

            snp_prevalence = ((alts>=alt_threshold)*(depths>D_min)).sum()
            snp_freq = snp_prevalence*1.0/(depths>D_min).sum()


            #if (population_freq==0) and (snp_freq==0):
            #    continue

            #if (population_freq==1) and (snp_freq==1):
            #    continue

            #if (snp_freq==0) or (snp_freq==1):
            #    continue


            depths_filter = depths[depths>D_min]

            freqs = alts[depths>D_min] / depths_filter
            samples_to_keep = samples[depths>D_min]


            f_max = max(freqs)
            # remove these sites because we don't examine them with the gamma
            if (f_max == 1) or (f_max == 0):
                continue

            for sample_i, freq_i in zip(samples_to_keep, freqs):

                pi_dict[location_aa][sample_i]['pi_sum_include_boundary'] += 2*freq_i*(1-freq_i)
                pi_dict[location_aa][sample_i]['n_sites_include_boundary'] += 1

                if (freq_i != 1) and (freq_i != 0):

                    pi_dict[location_aa][sample_i]['pi_sum_exclude_boundary'] += 2*freq_i*(1-freq_i)
                    pi_dict[location_aa][sample_i]['n_sites_exclude_boundary'] += 1



        for allowed_variant_type in allowed_variant_types:
            for sample in samples:

                if pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] > 0:
                    pi_dict[allowed_variant_type][sample]['pi_include_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_include_boundary']

                if pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] > 0:
                    pi_dict[allowed_variant_type][sample]['pi_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']


        intermediate_filename = intermediate_filename_template % species_name

        with open(intermediate_filename, 'wb') as handle:
            pickle.dump(pi_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def calculate_pi_from_snps_fmax_cutoff(f_max_range):

    intermediate_filename_template = config.data_directory+"pi_annotated_snps_f_max_cutoff/%s.dat"

    for species_name in good_species_list:

        sys.stderr.write("Loading whitelisted genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)

        snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        if os.path.isfile(snp_file_path) == False:
            continue

        # Open post-processed MIDAS output
        snp_file =  bz2.BZ2File(snp_file_path, "r")
        line = snp_file.readline() # header
        items = line.split()[1:]
        samples = numpy.array([item.strip() for item in items])


        pi_dict = {}

        for f_max_range_i in f_max_range:

            pi_dict[f_max_range_i] = {}

            for allowed_variant_type in allowed_variant_types:
                pi_dict[f_max_range_i][allowed_variant_type] = {}

                for sample in samples:
                    pi_dict[f_max_range_i][allowed_variant_type][sample] = {}

                    pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_sum_include_boundary'] = 0
                    pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_include_boundary'] = 0

                    pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_sum_exclude_boundary'] = 0
                    pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_exclude_boundary'] = 0


        sys.stderr.write("Iterating though SNPs...\n")
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

            # only examine core genes
            if gene_name not in core_genes:
                continue

            if location_aa not in allowed_variant_types:
                continue

            alts = []
            depths = []
            for item in items[1:]:
                subitems = item.split(",")
                alts.append(long(subitems[0]))
                depths.append(long(subitems[1]))
            alts = numpy.array(alts)
            depths = numpy.array(depths)
            refs = depths-alts


            # only try and predict prevalence if there's non-zero coverage for >= 20 hos
            if sum(depths>D_min) < min_n_samples:
                continue

            # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
            population_prevalence = ((alts>=refs)*(depths>D_min)).sum()
            population_freq = population_prevalence*1.0/(depths>D_min).sum()


            if population_freq>0.5:
                # alternate allele is in the majority
                # re-polarize for now
                alts,refs = refs,alts

            # Next calculate fraction of cohort where population minor allele is present at >=10% within-host frequency
            alt_threshold = numpy.ceil(depths*0.1)+0.5 #at least one read above 10%.

            snp_prevalence = ((alts>=alt_threshold)*(depths>D_min)).sum()
            snp_freq = snp_prevalence*1.0/(depths>D_min).sum()


            #if (population_freq==0) and (snp_freq==0):
            #    continue

            #if (population_freq==1) and (snp_freq==1):
            #    continue

            #if (snp_freq==0) or (snp_freq==1):
            #    continue


            depths_filter = depths[depths>D_min]

            freqs = alts[depths>D_min] / depths_filter
            samples_to_keep = samples[depths>D_min]


            f_max = max(freqs)
            # remove these sites because we don't examine them with the gamma
            if (f_max == 1) or (f_max == 0):
                continue

            for f_max_range_i in f_max_range:

                if f_max_range_i <= f_max:

                    for sample_i, freq_i in zip(samples_to_keep, freqs):

                        pi_dict[f_max_range_i][location_aa][sample_i]['pi_sum_include_boundary'] += 2*freq_i*(1-freq_i)
                        pi_dict[f_max_range_i][location_aa][sample_i]['n_sites_include_boundary'] += 1

                        if (freq_i != 1) and (freq_i != 0):

                            pi_dict[f_max_range_i][location_aa][sample_i]['pi_sum_exclude_boundary'] += 2*freq_i*(1-freq_i)
                            pi_dict[f_max_range_i][location_aa][sample_i]['n_sites_exclude_boundary'] += 1


        for f_max_range_i in f_max_range:
            for allowed_variant_type in allowed_variant_types:
                for sample in samples:

                    if pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_include_boundary'] > 0:
                        pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_include_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_include_boundary']

                    if pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_exclude_boundary'] > 0:
                        pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']


        intermediate_filename = intermediate_filename_template % species_name

        with open(intermediate_filename, 'wb') as handle:
            pickle.dump(pi_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_pi_dict(species_name):

    #if include_samples_with_strains == False:
    #    intermediate_filename_template = config.data_directory+"pi_annotated_snps_no_strains/%s.dat"

    #else:
    intermediate_filename_template = config.data_directory+"pi_annotated_snps/%s.dat"

    intermediate_filename = intermediate_filename_template % species_name

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b



def predict_prevalency():

    subject_sample_map = parse_HMP_data.parse_subject_sample_map()

    for species_name in good_species_list:

        pi_dict = load_pi_dict(species_name)

        sys.stderr.write("%s\n" % species_name)

        sys.stderr.write("Loading whitelisted genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)

        samples_in_dict_4D = [x for x in pi_dict['4D'].keys() if 'pi_exclude_boundary' in pi_dict['4D'][x] ]
        samples_in_dict_1D = [x for x in pi_dict['1D'].keys() if 'pi_exclude_boundary' in pi_dict['1D'][x] ]

        samples_in_dict_4D = numpy.asarray(samples_in_dict_4D)
        samples_in_dict_1D = numpy.asarray(samples_in_dict_1D)

        samples_from_dict = numpy.intersect1d(samples_in_dict_4D, samples_in_dict_1D)

        print(species_name)

        # Holds panel wide prevalence for each species
        #os.system('mkdir -p %ssnp_prevalences' % config.data_directory)
        snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        if os.path.isfile(snp_file_path) == False:
            continue


        # Open post-processed MIDAS output
        snp_file =  bz2.BZ2File(snp_file_path, "r")
        line = snp_file.readline() # header
        items = line.split()[1:]

        samples = numpy.array([item.strip() for item in items])
        samples_in_snp_and_pi_dict = numpy.intersect1d(samples, samples_from_dict)

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


        samples_idx = numpy.asarray([numpy.where(samples == sample)[0][0] for sample in samples_in_snp_and_pi_dict])
        samples_largest_clade_idx = numpy.asarray([numpy.where(samples == sample)[0][0] for sample in samples_in_snp_and_pi_dict if (sample in largest_clade_set)])

        samples_idx_dict = {}
        samples_idx_dict['all'] = samples_idx
        samples_idx_dict['largest_clade'] = samples_largest_clade_idx

        if len(samples_largest_clade_idx) < 20:
            continue

        things_to_measure = ['f_max', 'observed_prevalence', 'predicted_prevalence', 'f_mean', 'f_var', 'f_mean_no_zeros', 'f_var_no_zeros', 'observed_prevalence/predicted_prevalence', 'predicted_f_mean']

        # fix this mess.......
        pi_samples_dict = {}
        #count_variant_type_dict = {}
        #site_occupancy_dict = {}
        freq_dict = {}
        for clade_type in clade_types:
            freq_dict[clade_type] = {}
            pi_samples_dict[clade_type] = {}

            for allowed_variant_type in allowed_variant_types:
                freq_dict[clade_type][allowed_variant_type] = {}
                for k in things_to_measure:
                    freq_dict[clade_type][allowed_variant_type][k] = []

                pi_samples = [pi_dict[allowed_variant_type][sample]['pi_exclude_boundary'] for sample in samples[samples_idx_dict[clade_type]]]
                pi_samples = numpy.asarray(pi_samples)
                pi_samples_dict[clade_type][allowed_variant_type] = pi_samples

        # add samples to the dict samples_idx_dict
        freq_dict['all']['samples'] = samples_in_snp_and_pi_dict.tolist()
        freq_dict['largest_clade']['samples'] = [sample for sample in samples_in_snp_and_pi_dict if (sample in largest_clade_set)]

        sys.stderr.write("Calculating SNP prevalences...\n")
        num_sites_processed = 0
        for line in snp_file:

            num_sites_processed+=1
            if num_sites_processed%50000==0:
                sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))
                #if debug:
                #    break

            items = line.split()
            # Load information about site
            info_items = items[0].split("|")
            chromosome = info_items[0]
            location = long(info_items[1])
            gene_name = info_items[2]
            location_aa = info_items[3]

            # only examine core genes
            if gene_name not in core_genes:
                continue

            if location_aa not in allowed_variant_types:
                continue

            #if (count_variant_type_dict['1D'] >= count) and (count_variant_type_dict['4D'] >= count):
            #    break

            #while (count_variant_type_dict['1D'] < count) and (count_variant_type_dict['4D'] < count):
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


            for clade_type in clade_types:

                if clade_type == 'largest_clade':
                    alts_iter = alts[samples_largest_clade_idx]
                    depths_iter = depths[samples_largest_clade_idx]
                    refs_iter = refs[samples_largest_clade_idx]

                else:
                    alts_iter = alts[samples_idx]
                    depths_iter = depths[samples_idx]
                    refs_iter = refs[samples_idx]


                # only look at alleles where the coveragae threshold is met in 20 hosts
                if sum(depths>D_min) < min_n_samples:
                    continue

                # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
                population_prevalence = ((alts_iter>=refs_iter)*(depths_iter>D_min)).sum()
                population_freq = population_prevalence*1.0/(depths_iter>D_min).sum()

                if population_freq>0.5:
                    # alternate allele is in the majority
                    # re-polarize for now
                    alts_iter,refs_iter = refs_iter,alts_iter

                # Next calculate fraction of cohort where population minor allele is present at >=10% within-host frequency
                alt_threshold = numpy.ceil(depths_iter*0.1)+0.5 #at least one read above 10%.

                snp_prevalence = ((alts_iter>=alt_threshold)*(depths_iter>D_min)).sum()
                snp_freq = snp_prevalence*1.0/(depths_iter>D_min).sum()

                depths_filter = depths_iter[depths_iter>D_min]

                pi_samples_filter = pi_samples_dict[clade_type][location_aa][depths_iter>D_min]

                freqs = alts_iter[depths_iter>D_min] / depths_filter

                # make sure you have enough samples after filtering on coverage
                if len(freqs) < min_n_samples:
                    continue

                f_max = max(freqs)

                # absorbing boundary, ignore
                if (f_max == 1) or (f_max == 0):
                    continue

                prevalence = len(freqs[freqs>0]) / len(freqs)

                if prevalence == 0:
                    continue

                freq_dict[clade_type][location_aa]['f_max'].append(f_max)
                freq_dict[clade_type][location_aa]['observed_prevalence'].append(prevalence)

                freq_dict[clade_type][location_aa]['f_mean'].append(numpy.mean(freqs))
                freq_dict[clade_type][location_aa]['f_var'].append(numpy.var(freqs))

                freq_dict[clade_type][location_aa]['f_mean_no_zeros'].append(numpy.mean(freqs[freqs>0]))
                freq_dict[clade_type][location_aa]['f_var_no_zeros'].append(numpy.var(freqs[freqs>0]))


                predicted_prevalence = 1 - numpy.mean( (1 + f_max * depths_filter) ** (-1*pi_samples_filter) )

                freq_dict[clade_type][location_aa]['predicted_prevalence'].append(predicted_prevalence)
                freq_dict[clade_type][location_aa]['observed_prevalence/predicted_prevalence'].append(prevalence/predicted_prevalence)

                freq_dict[clade_type][location_aa]['predicted_f_mean'].append( f_max * numpy.mean(pi_samples_filter) )

                #count_variant_type_dict[location_aa] += 1

        snp_file.close()
        sys.stderr.write("Done!\n")

        intermediate_filename_template = config.data_directory+"predicted_observed_prevalence/%s_%s.dat"

        for clade_type in clade_types:

            intermediate_filename = intermediate_filename_template % (species_name, clade_type)

            with open(intermediate_filename, 'wb') as handle:
                pickle.dump(freq_dict[clade_type], handle, protocol=pickle.HIGHEST_PROTOCOL)




def subsample_data(n_subsample=1000, n_iterations=10000, n_to_plot=10000):

    subsample_dict = {}

    for species_name in good_species_list:

        for clade_type in clade_types:

            intermediate_filename = '%spredicted_observed_prevalence/%s_%s.dat'% (config.data_directory, species_name, clade_type)

            if os.path.isfile(intermediate_filename) == False:
                continue

            with open(intermediate_filename, 'rb') as handle:
                snp_dict = pickle.load(handle)

            if species_name not in subsample_dict:

                subsample_dict[species_name] = {}

            subsample_dict[species_name][clade_type] = {}

            subsample_dict[species_name][clade_type]['samples'] = snp_dict['samples']

            for allowed_variant_type in allowed_variant_types:

                predicted_prevalence = snp_dict[allowed_variant_type]['predicted_prevalence']
                predicted_prevalence = numpy.asarray(predicted_prevalence)

                observed_prevalence = snp_dict[allowed_variant_type]['observed_prevalence']
                observed_prevalence = numpy.asarray(observed_prevalence)

                prevalence_ratio = snp_dict[allowed_variant_type]['observed_prevalence/predicted_prevalence']
                prevalence_ratio = numpy.asarray(prevalence_ratio)

                f_mean = snp_dict[allowed_variant_type]['f_mean']
                f_mean = numpy.asarray(f_mean)
                f_mean_log10 = numpy.log10(f_mean)

                predicted_f_mean = snp_dict[allowed_variant_type]['predicted_f_mean']
                predicted_f_mean = numpy.asarray(predicted_f_mean)

                f_max = snp_dict[allowed_variant_type]['f_max']
                f_max = numpy.asarray(f_max)
                f_max_log10 = numpy.log10(f_max)

                if len(predicted_prevalence) < n_to_plot*2:
                    continue

                # make dictionary entry
                subsample_dict[species_name][clade_type][allowed_variant_type] = {}

                idx_to_plot = numpy.random.choice(len(predicted_prevalence), size=n_to_plot, replace=False)

                subsample_dict[species_name][clade_type][allowed_variant_type]['predicted_prevalence_to_plot'] = predicted_prevalence[idx_to_plot].tolist()
                subsample_dict[species_name][clade_type][allowed_variant_type]['observed_prevalence_to_plot'] = observed_prevalence[idx_to_plot].tolist()

                subsample_dict[species_name][clade_type][allowed_variant_type]['f_mean_to_plot'] = f_mean[idx_to_plot].tolist()
                subsample_dict[species_name][clade_type][allowed_variant_type]['predicted_f_mean_to_plot'] = predicted_f_mean[idx_to_plot].tolist()

                subsample_dict[species_name][clade_type][allowed_variant_type]['mean_prevalence_ratio'] = numpy.mean(prevalence_ratio)
                abs_error = numpy.absolute(predicted_prevalence - observed_prevalence)
                relative_error = abs_error / observed_prevalence
                mre = numpy.mean(relative_error)
                subsample_dict[species_name][clade_type][allowed_variant_type]['relative_error_to_plot'] = relative_error[idx_to_plot].tolist()

                subsample_dict[species_name][clade_type][allowed_variant_type]['f_max_to_plot'] = f_max[idx_to_plot].tolist()


                mean_error = numpy.mean(predicted_prevalence - observed_prevalence)

                mae = numpy.mean(abs_error)
                mse = numpy.mean((predicted_prevalence - observed_prevalence)**2)
                std_error_ae = numpy.std(abs_error) / numpy.sqrt(len(abs_error))


                # f mean
                abs_error_f_mean = numpy.absolute(predicted_f_mean - f_mean)
                relative_error_f_mean = abs_error_f_mean / f_mean
                mre_f_mean = numpy.mean(relative_error_f_mean)


                def bootstrap_sample(error_):

                    boostrap_ = []

                    for i in range(n_iterations):

                        error_i = numpy.random.choice(error_, size=n_subsample, replace=True)
                        mean_error_i = numpy.mean(error_i)
                        boostrap_.append(mean_error_i)

                    boostrap_ = numpy.asarray(boostrap_)
                    boostrap_ = numpy.sort(boostrap_)

                    boostrap_025 = boostrap_[int(n_iterations*0.025)]
                    boostrap_975 = boostrap_[int(n_iterations*0.975)]

                    return boostrap_025, boostrap_975


                mae_boostrap_025, mae_boostrap_975 = bootstrap_sample(abs_error)
                mre_boostrap_025, mre_boostrap_975 = bootstrap_sample(relative_error)


                subsample_dict[species_name][clade_type][allowed_variant_type]['ME'] = mean_error
                subsample_dict[species_name][clade_type][allowed_variant_type]['AE_SE'] = std_error_ae

                subsample_dict[species_name][clade_type][allowed_variant_type]['MAE'] = mae
                subsample_dict[species_name][clade_type][allowed_variant_type]['MAE_025'] = mae_boostrap_025
                subsample_dict[species_name][clade_type][allowed_variant_type]['MAE_975'] = mae_boostrap_975

                subsample_dict[species_name][clade_type][allowed_variant_type]['MSE'] = mse

                subsample_dict[species_name][clade_type][allowed_variant_type]['MRE'] = mre
                subsample_dict[species_name][clade_type][allowed_variant_type]['MRE_025'] = mre_boostrap_025
                subsample_dict[species_name][clade_type][allowed_variant_type]['MRE_975'] = mre_boostrap_975


                subsample_dict[species_name][clade_type][allowed_variant_type]['MRE_f_mean'] = mre_f_mean

                # relationship between error and f_max
                hist, bin_edges = numpy.histogram(f_max_log10, density=True, bins=40)

                bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]

                prevalence_predicted_mean = []
                relative_error_mean = []
                for i in range(0, len(bin_edges)-1):
                    idx_i = (f_max_log10 > bin_edges[i]) & (f_max_log10 < bin_edges[i+1])
                    prevalence_predicted_mean.append(numpy.mean(numpy.log10(predicted_prevalence[idx_i])))
                    relative_error_mean.append(numpy.mean(relative_error[idx_i]))

                subsample_dict[species_name][clade_type][allowed_variant_type]['f_max_line'] = bins_mean
                subsample_dict[species_name][clade_type][allowed_variant_type]['f_max_vs_predicted_prevalence_line'] = prevalence_predicted_mean
                subsample_dict[species_name][clade_type][allowed_variant_type]['f_max_vs_relative_error_line'] = relative_error_mean


                # relationship between error and mean f
                relative_error_mean = []
                hist, bin_edges = numpy.histogram(f_mean_log10, density=True, bins=40)
                bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]
                for i in range(0, len(bin_edges)-1):
                    idx_i = (f_mean_log10 > bin_edges[i]) & (f_mean_log10 < bin_edges[i+1])
                    relative_error_mean.append(numpy.mean(relative_error[idx_i]))

                subsample_dict[species_name][clade_type][allowed_variant_type]['f_mean_line'] = bins_mean
                subsample_dict[species_name][clade_type][allowed_variant_type]['f_mean_vs_relative_error_line'] = relative_error_mean


                estimates = ['f_max', 'f_mean']
                # regress error against these estimates

                for e in estimates:

                    e_ = snp_dict[allowed_variant_type][e]
                    e_ = numpy.asarray(e_)

                    slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(e_), numpy.log10(relative_error))

                    subsample_dict[species_name][clade_type][allowed_variant_type][e] = {}

                    subsample_dict[species_name][clade_type][allowed_variant_type][e]['to_plot'] =  e_[idx_to_plot].tolist()
                    subsample_dict[species_name][clade_type][allowed_variant_type][e]['slope'] = slope
                    subsample_dict[species_name][clade_type][allowed_variant_type][e]['intercept'] = intercept
                    subsample_dict[species_name][clade_type][allowed_variant_type][e]['r_value'] = r_value
                    subsample_dict[species_name][clade_type][allowed_variant_type][e]['p_value'] = p_value
                    subsample_dict[species_name][clade_type][allowed_variant_type][e]['std_err'] = std_err

                    subsample_dict[species_name][clade_type][allowed_variant_type][e]['mean'] = numpy.mean(e_)


    intermediate_filename = config.data_directory+"predicted_prevalence.dat"

    with open(intermediate_filename, 'wb') as handle:
        pickle.dump(subsample_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_predicted_prevalence_subsample_dict():

    intermediate_filename = config.data_directory+"predicted_prevalence.dat"

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b



if __name__=='__main__':

    #calculate_pi_from_snps_fmax_cutoff(f_max_range)

    #calculate_pi_from_snps()
    #predict_prevalency()
    subsample_data()
