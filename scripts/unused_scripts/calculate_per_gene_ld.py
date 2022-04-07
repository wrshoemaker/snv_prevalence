#import matplotlib
#matplotlib.use('Agg')
import parse_midas_data
#import pylab
#from pylab import *
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
import os
#import pandas
import parse_patric
import pickle
import sample_utils
import parse_HMP_data

import calculate_linkage_disequilibria

# code from plot_kegg_pi_distribution.py

variant_types = ['4D','1D']

bin_width_exponent = 1.1


################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize



sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

#good_species_list = parse_midas_data.parse_good_species_list()
good_species_list = ['Bacteroides_vulgatus_57955']
ld_dict_all_species = calculate_linkage_disequilibria.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)
#print(ld_dict_all_species.keys())

ld_species = list(ld_dict_all_species.keys())

print(ld_dict_all_species['Bacteroides_vulgatus_57955']['4D'].keys())
