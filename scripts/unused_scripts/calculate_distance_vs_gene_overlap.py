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

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
from math import log10,ceil,fabs
from numpy.random import randint, choice

from itertools import combinations

import midas_db_utils
import parse_HMP_data


min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates

# code from plot_closely_related_strains.py


min_sample_size = config.between_host_min_sample_size

intermediate_filename_template = '%s%s.txt.gz'

pairwise_directory = '%spairwise_distance/' % (parse_midas_data.data_directory)

################################################################################
#
# Standard header to read in argument information
#
################################################################################
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
sys.stderr.write("Done!\n")

good_species_list = parse_midas_data.parse_good_species_list()
if species!='all':
    good_species_list = [species]
if debug and len(good_species_list)>3.5:
    #good_species_list = good_species_list[:3]
    good_species_list = ['Bacteroides_vulgatus_57955']


for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
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


    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough hosts!\n")
        continue
    else:
        sys.stderr.write("Found %d unique hosts!\n" % len(snp_samples))


    # Load divergence matrices
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, '4D', allowed_samples=snp_samples)
    snp_samples = numpy.array(dummy_samples)
    substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    sys.stderr.write("Done!\n")


    # get gene coverage

    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
    sys.stderr.write("Done!\n")

    sys.stderr.write("Loaded gene info for %d samples\n" % len(gene_samples))


    # Calculate matrix of number of genes that differ
    sys.stderr.write("Calculating matrix of gene differences...\n")
    gene_difference_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_reads_matrix, gene_depth_matrix, marker_coverages)

    # Now need to make the gene samples and snp samples match up
    desired_samples = gene_samples[marker_coverages>min_coverage]



    # Calculate which pairs of idxs belong to the same sample, which to the same subject
    # and which to different subjects
    desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = sample_utils.calculate_subject_pairs(subject_sample_map, desired_samples)

    snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)

    same_sample_snp_idxs  = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  desired_same_sample_idxs)
    same_sample_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_sample_idxs)

    same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)
    same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)

    diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)
    diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)

    #typical_same_subject_snp_opportunities = numpy.median(snp_opportunity_matrix[same_subject_snp_idxs])
    #typical_diff_subject_snp_opportunities = numpy.median(snp_opportunity_matrix[diff_subject_snp_idxs])

    #Lsnps = typical_diff_subject_snp_opportunities
    #Lgenes = typical_diff_subject_gene_opportunities


    diff_subject_snp_plowers = []
    diff_subject_snp_puppers = []
    diff_subject_gene_plowers = []
    diff_subject_gene_puppers = []

    diff_subject_snp_expected = []
    diff_subject_gene_expected = []

    record_strs = [", ".join(['Sample_1_idx', 'Sample_2_idx', 'snp_pexpected', 'snp_plower', 'snp_pupper', 'gene_pexpected', 'gene_plower', 'gene_pupper'])]

    for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):

        i = diff_subject_snp_idxs[0][sample_pair_idx]
        j = diff_subject_snp_idxs[1][sample_pair_idx]

        #sample_name_1 = snp_samples[i]
        #sample_name_2 = snp_samples[j]

        snp_plower,snp_pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[i,j], snp_opportunity_matrix[i,j])

        #diff_subject_snp_plowers.append(plower)
        #diff_subject_snp_puppers.append(pupper)

        snp_pexpected = snp_difference_matrix[i,j] / snp_opportunity_matrix[i,j]
        #diff_subject_snp_expected.append(snp_difference_matrix[i,j] / snp_opportunity_matrix[i,j])



        i = diff_subject_gene_idxs[0][sample_pair_idx]
        j = diff_subject_gene_idxs[1][sample_pair_idx]
        gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages)

        gene_plower,gene_pupper = stats_utils.calculate_poisson_rate_interval(gene_difference_matrix[i,j], gene_opportunity_matrix[i,j],alpha)

        gene_pexpected = gene_difference_matrix[i,j] / gene_opportunity_matrix[i,j]

        #diff_subject_gene_plowers.append(plower)
        #diff_subject_gene_puppers.append(pupper)

        #diff_subject_gene_expected.append(gene_difference_matrix[i,j] / gene_opportunity_matrix[i,j])

        record_str = ", ".join([str(i), str(j), str(snp_pexpected), str(snp_plower), str(snp_pupper), str(gene_pexpected), str(gene_plower), str(gene_pupper)])
        record_strs.append(record_str)


    sys.stderr.write("Done with %s!\n" % species_name)

    sys.stderr.write("Writing intermediate file...\n")
    intermediate_filename = intermediate_filename_template % (pairwise_directory, species_name)
    file = gzip.open(intermediate_filename,"w")
    record_str = "\n".join(record_strs)
    file.write(record_str)
    file.close()
    sys.stderr.write("Done!\n")



    #snp_sample_pairs_idx = combinations(list(range(len(snp_samples))),2)

    #record_strs = [", ".join(['Sample_1', 'Sample_2', 'Substitiot rate', 'Gene overlap'])]

    #for snp_sample_pair_idx in snp_sample_pairs_idx:

    #    snp_sample_1_idx = snp_sample_pair_idx[0]
    #    snp_sample_2_idx = snp_sample_pair_idx[1]

    #    sample_1 = snp_samples[snp_sample_1_idx]
    #    sample_2 = snp_samples[snp_sample_2_idx]

    #    genetic_distance = substitution_rate[snp_sample_1_idx, snp_sample_2_idx]

    #    genes_sum = gene_presence_matrix[:,snp_sample_1_idx] + gene_presence_matrix[:,snp_sample_2_idx]
    #    genes_intersect = genes_sum[genes_sum ==2]

    #    gene_interect_proportion = len(genes_intersect) / len(genes_sum)

    #    record_str = ", ".join([sample_1, sample_2, str(genetic_distance), str(gene_interect_proportion)])

    #    record_strs.append(record_str)


    #sys.stderr.write("Done with %s!\n" % species_name)

    #sys.stderr.write("Writing intermediate file...\n")
    #intermediate_filename = intermediate_filename_template % (pairwise_directory, species_name)
    #file = gzip.open(intermediate_filename,"w")
    #record_str = "\n".join(record_strs)
    #file.write(record_str)
    #file.close()
    #sys.stderr.write("Done!\n")
