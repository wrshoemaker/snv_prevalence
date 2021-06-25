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

good_species_list = parse_midas_data.parse_good_species_list()


min_sample_size = config.between_host_min_sample_size



subject_sample_map = parse_HMP_data.parse_subject_sample_map()

for species_name in good_species_list:


    sys.stderr.write("Loading whitelisted genes...\n")
    core_genes = core_gene_utils.parse_core_genes(species_name)


    print(species_name)


    # get samples from largest clade
    snp_samples = diversity_utils.calculate_haploid_samples(species_name)
    if len(snp_samples) < min_sample_size:
        continue
    snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
    # get clade idxs
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
    #snp_samples = numpy.array(dummy_samples)
