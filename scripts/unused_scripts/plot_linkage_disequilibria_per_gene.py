import sys, os
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria
import calculate_linkage_disequilibria_per_gene

import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#import calculate_snv_distances
import figure_utils
from math import log10,ceil
#import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.random import randint, multinomial
import parse_HMP_data

import figure_utils
import ld_utils

import matplotlib as mpl


variant_types = ['4D','1D']

bin_width_exponent = 1.1



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

ld_dict_all_species = calculate_linkage_disequilibria_per_gene.calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=bin_width_exponent)

ld_species = list(ld_dict_all_species.keys())

n_cols = 6
species_nested_list  = list(ld_utils.chunks(ld_species, n_cols))
n_rows = len(species_nested_list)


#for species in ld_species:
