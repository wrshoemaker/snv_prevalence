from __future__ import division
import matplotlib
matplotlib.use('Agg')
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import gene_diversity_utils
import stats_utils
import parse_HMP_data

import pickle
################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size

# Minimum median coverage of sample to look at
min_coverage = 20
min_passed_sites_per_gene=10
min_passed_sites_per_person=100


################################################################################
#Bacteroides_vulgatus_57955

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

# Load gene presence/absence information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Loaded %d genes across %d samples\n" % gene_depth_matrix.shape)
sys.stderr.write("Done!\n")

min_marker_coverage = 20
high_coverage_samples = samples[marker_coverages>=min_marker_coverage]
sys.stderr.write("Focusing on %d high coverage samples...\n" % len(high_coverage_samples))

# Load metaphlan2 genes
#sys.stderr.write("Loading metaphlan2 genes...\n")
#metaphlan2_genes = set(parse_midas_data.load_metaphlan2_genes(species_name))
#metaphlan2_gene_idxs = numpy.array([gene_name in metaphlan2_genes for gene_name in gene_names])
#sys.stderr.write("Done! (%d genes)\n" % len(metaphlan2_genes))

# Load reference genes
sys.stderr.write("Loading reference genes...\n")
reference_genes = set(parse_midas_data.load_reference_genes(species_name))
reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])
sys.stderr.write("Done! (%d genes)\n" % len(reference_genes))

#core genes
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! %d core genes\n" % len(core_genes))




######################
# Load coverage data #
######################

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

###############################################################
# Compute Pi within patients to figure out which are haploid  #
###############################################################

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
#samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_variant_types=set(['4D']), allowed_genes=core_genes, debug=debug)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi_new(species_name, allowed_variant_types=set(['4D']), allowed_genes=core_genes, debug=debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities


######################
# compute median cov #
######################
median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

###############################################################
# Indexes for SNP samples that have high coverage             #
###############################################################

# Only plot samples above a certain depth threshold that are "haploids"
low_pi_snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
high_pi_snp_samples = samples[(median_coverages>=min_coverage)*(pis>1e-03)]


##########################################################
# load SNP info
##########################################################
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug, allowed_samples=low_pi_snp_samples)
sys.stderr.write("Done!\n")

#########################################################
# figure out which samples belong to the same subject   #
#########################################################
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)









gene_copy_numbers = (gene_depth_matrix/marker_coverages)[:,marker_coverages>=min_marker_coverage]
gene_presence_calls = gene_copy_numbers>0.1

gene_prevalence_map = {gene_names[i]: gene_presence_calls[i,:].sum() for i in xrange(0,len(gene_names))}

gene_prevalences = (gene_copy_numbers>0.1).sum(axis=1) / gene_copy_numbers.shape[1]

bin_integer = 0.10

prevalence_bin_dict = {}
prevalence_bins = numpy.linspace(0, 1, num=11)

for i in range(len(prevalence_bins)-1):

    #prevalence_bin_dict[(prevalence_bins[i], prevalence_bins[i+1])] = {}
    prevalence_bin_genes = gene_names[(gene_prevalences >= prevalence_bins[i]) & (gene_prevalences < prevalence_bins[i+1])]
    #prevalence_bin_dict[(prevalence_bins[i], prevalence_bins[i+1])]['gene_name'] = prevalence_bin_genes

    prevalence_bin_genes_set = set(prevalence_bin_genes)

    allowed_genes_prevalence = prevalence_bin_genes_set & set(passed_sites_map.keys())


    ##########################################################
    # compute total fixations, genome-wide core genes
    ##########################################################
    # Calculate fixation matrices
    min_change=0.8
    sys.stderr.write("Calculating 4D fixation matrix...\n")
    fixation_matrix_syn, fixation_opportunities_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']),allowed_genes=allowed_genes_prevalence, min_change=min_change)
    sys.stderr.write("Calculating 1D fixation matrix...\n")
    fixation_matrix_non, fixation_opportunities_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']),allowed_genes=allowed_genes_prevalence, min_change=min_change)
    sys.stderr.write("Calculating total fixation matrix...\n")
    fixation_matrix_all, fixation_opportunities_all = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_genes=allowed_genes_prevalence, min_change=min_change)

    sys.stderr.write("Done!\n")

    # Calculate fraction nonsynonymous
    dN = fixation_matrix_non/fixation_opportunities_non
    dS = fixation_matrix_syn/fixation_opportunities_syn
    dNplusdS = (dN+dS)
    fraction_nonsynonymous = dN/(dNplusdS+(dNplusdS==0))
    fraction_nonsynonymous_flatten = fraction_nonsynonymous[numpy.triu_indices(fraction_nonsynonymous.shape[1], 1)]
    fraction_nonsynonymous_flatten = fraction_nonsynonymous_flatten[((~numpy.isnan(fraction_nonsynonymous_flatten)) | (fraction_nonsynonymous_flatten >0))]

    #print(len(fraction_nonsynonymous_flatten))

    #fraction_nonsynonymous_flatten = numpy.triu(fraction_nonsynonymous, 1)

    #fraction_nonsynonymous[numpy.triu_indices(fraction_nonsynonymous.shape[1], 1)]
    #fraction_nonsynonymous_flatten = fraction_nonsynonymous_flatten[((~numpy.isnan(fraction_nonsynonymous_flatten)) | (fraction_nonsynonymous_flatten >0))]

    prevalence_bin_dict[(prevalence_bins[i], prevalence_bins[i+1])] = {}
    #prevalence_bin_dict[(prevalence_bins[i], prevalence_bins[i+1])][] = numpy.mean(fraction_nonsynonymous_flatten)


    pi_matrix_syn, avg_pi_matrix_syn, passed_sites_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=allowed_genes_prevalence)

    avg_pi_matrix_syn = avg_pi_matrix_syn/(passed_sites_syn+(passed_sites_syn==0))
    avg_pi_matrix_syn = numpy.clip(avg_pi_matrix_syn,1e-06,1)

    avg_pi_matrix_flatten = avg_pi_matrix_syn[numpy.triu_indices(avg_pi_matrix_syn.shape[1], 1)]
    avg_pi_matrix_avg = numpy.mean(avg_pi_matrix_flatten)
    avg_pi_matrix_se = numpy.std(avg_pi_matrix_flatten) / numpy.sqrt(len(avg_pi_matrix_flatten))

    prevalence_bin_dict[(prevalence_bins[i], prevalence_bins[i+1])]['avg_pi_matrix_avg'] = avg_pi_matrix_avg
    prevalence_bin_dict[(prevalence_bins[i], prevalence_bins[i+1])]['avg_pi_matrix_se'] = avg_pi_matrix_se




saved_data_file='%sprevalence_pi/prevalence_pi_%s.dat' % (parse_midas_data.data_directory, species_name)
with open(saved_data_file, 'wb') as outfile:
    pickle.dump(prevalence_bin_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)

#stringent_gene_prevalences = (gene_copy_numbers>0.5).sum(axis=1)
#lax_gene_prevalences = (gene_copy_numbers>0.01).sum(axis=1)

#prevalence_bins = numpy.arange(0,gene_copy_numbers.shape[1]+2)-0.5
#prevalences = prevalence_bins[1:]+(prevalence_bins[1]-prevalence_bins[0])/2
#prevalences /= prevalences[-1]

#prevalence_counts = numpy.histogram(gene_prevalences, prevalence_bins)[0]*1.0
#prevalence_counts /= prevalence_counts.sum()
