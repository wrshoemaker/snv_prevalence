

import config
import parse_midas_data

import os
import parse_HMP_data

import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
import sample_utils

##################
# print two files
# 1) fraction of snps in each bin for all 4D and 1D sites for each species, all in one file
# 2) bin identity and AA identity for each site, do this for *each* species
##################


min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_divergence_threshold = 5e-04
min_change = 0.8
allowed_variant_types = set(['1D','2D','3D','4D'])
maf_cutoff = 0.5
#n = 30
min_sample_size = config.between_host_min_sample_size


#def get_maf(species_name):


#sfs_file_path = '%smaf_all_species.txt' % (parse_midas_data.data_directory)


# Oscillibacter_sp_60799 error ?/?




if __name__=='__main__':


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("--species", help="Name of specific species to run code on", default="all")

    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species=args.species


    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sample_country_map = parse_HMP_data.parse_sample_country_map()
    sys.stderr.write("Done!\n")

    good_species_list = parse_midas_data.parse_good_species_list()
    if species!='all':
        good_species_list = [species]
    if debug and len(good_species_list)>3.5:
        #good_species_list = good_species_list[:1]
        #good_species_list = ['Bacteroides_uniformis_57318']
        good_species_list = ['Bacteroides_vulgatus_57955']



    # header of the output file.
    record_strs = [", ".join(['Species', 'SNPSamples', 'MAFbins', 'MAF1D', 'MAF4D'])]


    for species_name in good_species_list:

        # Only plot samples above a certain depth threshold that are "haploids"
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        else:
            sys.stderr.write("Found %d haploid samples!\n" % len(snp_samples))


        sys.stderr.write("Calculating unique hosts...\n")
        # Only consider one sample per person
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        snp_samples = [s.decode("utf-8")  for s in snp_samples]

        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough hosts!\n")
            continue
        else:
            sys.stderr.write("Found %d unique hosts!\n" % len(snp_samples))


        ###############
        #
        # Do calculation!
        #
        ###############

        sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        sys.stderr.write("Calculating matrix...\n")
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = numpy.array(dummy_samples)
        substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        sys.stderr.write("Done!\n")



        #################
        #
        # Cluster samples into clades based on distance matrix
        #
        #################
        sys.stderr.write("Clustering samples with low divergence...\n")
        coarse_grained_idxs, coarse_grained_cluster_list = clade_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold)

        coarse_grained_samples = snp_samples[coarse_grained_idxs]
        clade_sets = clade_utils.load_manual_clades(species_name)

        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)

        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])

        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]


        sys.stderr.write("Top level: %d clades, %s\n" % (len(clade_sets), str(clade_sizes)))
        sys.stderr.write("Max: %d\n" % len(largest_clade_samples))


        if len(largest_clade_samples)<3:
            sys.stderr.write("Too few samples! Quitting...\n")
            sys.exit(1)

        sys.stderr.write("Continuing with %d samples...\n" % len(largest_clade_samples))

        # Load SNP information for species_name
        sys.stderr.write("Re-loading %s...\n" % species_name)

        sys.stderr.write("Loading core genes...\n")
        core_genes = parse_midas_data.load_core_genes(species_name)
        sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))


        snp_difference_matrix = numpy.array([])
        snp_opportunity_matrix = numpy.array([])

        synonymous_difference_matrix = numpy.array([])
        synonymous_opportunity_matrix = numpy.array([])

        nonsynonymous_difference_matrix = numpy.array([])
        nonsynonymous_opportunity_matrix = numpy.array([])


        maf_bins = []
        mafs = []

        count_bins = []
        count_locations = []

        synonymous_sfs = []
        nonsynonymous_sfs = []

        synonymous_count_sfs = []
        nonsynonymous_count_sfs = []

        synonymous_pi_weighted_counts = 0
        nonsynonymous_pi_weighted_counts = 0


        final_line_number = 0
        while final_line_number >= 0:

            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_samples=largest_clade_samples,allowed_genes=core_genes, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))

            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of snp differences...\n")
            chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, min_change=min_change)
            sys.stderr.write("Done!\n")



            if snp_difference_matrix.shape[0]==0:
                snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
                snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                synonymous_difference_matrix = numpy.zeros_like(snp_difference_matrix)
                synonymous_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)
                nonsynonymous_difference_matrix = numpy.zeros_like(snp_difference_matrix)
                nonsynonymous_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)

                n = len(snp_samples)

                maf_bins = numpy.arange(1,n+1)*1.0/n
                maf_bins -= (maf_bins[1]-maf_bins[0])/2
                maf_bins[0]=-0.1
                maf_bins[-1] = 1.1
                mafs = numpy.arange(1,n)*1.0/n

                count_bins = numpy.arange(1,n+1)-0.5
                count_locations = numpy.arange(1,n)

                synonymous_sfs = numpy.zeros_like(mafs)
                nonsynonymous_sfs = numpy.zeros_like(mafs)

                synonymous_count_sfs = numpy.zeros_like(count_locations)
                nonsynonymous_count_sfs = numpy.zeros_like(count_locations)


            snp_difference_matrix += chunk_snp_difference_matrix
            snp_opportunity_matrix += chunk_snp_opportunity_matrix

            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of 4D differences...\n")
            chunk_synonymous_difference_matrix, chunk_synonymous_opportunity_matrix =    diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, min_change=min_change,allowed_variant_types=set(['4D']))
            sys.stderr.write("Done!\n")

            synonymous_difference_matrix += chunk_synonymous_difference_matrix
            synonymous_opportunity_matrix += chunk_synonymous_opportunity_matrix

            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of 1D differences...\n")
            chunk_nonsynonymous_difference_matrix, chunk_nonsynonymous_opportunity_matrix =    diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, min_change=min_change,allowed_variant_types=set(['1D']))
            sys.stderr.write("Done!\n")

            nonsynonymous_difference_matrix += chunk_nonsynonymous_difference_matrix
            nonsynonymous_opportunity_matrix += chunk_nonsynonymous_opportunity_matrix

            sys.stderr.write("Calculating the SFS...\n")
            #chunk_synonymous_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_variant_types = set(['4D']), allowed_genes=core_genes)
            #chunk_nonsynonymous_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_variant_types = set(['1D']), allowed_genes=core_genes)

            #chunk_synonymous_sfs, dummy = numpy.histogram(chunk_synonymous_freqs, bins=maf_bins)
            #synonymous_sfs += chunk_synonymous_sfs

            #chunk_nonsynonymous_sfs, dummy = numpy.histogram(chunk_nonsynonymous_freqs, bins=maf_bins)
            #nonsynonymous_sfs += chunk_nonsynonymous_sfs

            sys.stderr.write("Calculating count SFS...\n")
            #chunk_synonymous_counts, chunk_synonymous_weights = diversity_utils.calculate_pooled_counts(allele_counts_map, passed_sites_map, allowed_variant_types = set(['4D']), allowed_genes=core_genes,pi_min_k=4)
            #chunk_nonsynonymous_counts, chunk_nonsynonymous_weights = diversity_utils.calculate_pooled_counts(allele_counts_map, passed_sites_map, allowed_variant_types = set(['1D']), allowed_genes=core_genes,pi_min_k=4)

            #chunk_synonymous_count_sfs, dummy = numpy.histogram(chunk_synonymous_counts, bins=count_bins)
            #synonymous_count_sfs += chunk_synonymous_count_sfs

            #chunk_nonsynonymous_count_sfs, dummy = numpy.histogram(chunk_nonsynonymous_counts, bins=count_bins)
            #nonsynonymous_count_sfs += chunk_nonsynonymous_count_sfs

            #synonymous_pi_weighted_counts += chunk_synonymous_weights
            #nonsynonymous_pi_weighted_counts += chunk_nonsynonymous_weights

            diversity_utils.calculate_prevalence_matrix(allele_counts_map, passed_sites_map, allowed_variant_types = set(['1D']), allowed_genes=core_genes,pi_min_k=4)
