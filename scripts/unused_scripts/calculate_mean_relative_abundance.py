from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path
import random
from collections import Counter

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import plot_utils
import prevalence_utils

from scipy.stats import gamma, gaussian_kde
import scipy.stats as stats

import calculate_predicted_prevalence_mapgd

species_list = prevalence_utils.species_to_run


samples_dict = {}
for species_name in species_list:

    pi_dict_species = calculate_predicted_prevalence_mapgd.load_pi_dict(species_name)
    samples = list(pi_dict_species['4D'].keys())
    samples_dict[species_name] = samples


depth_file = bz2.BZ2File("%sspecies/relative_abundance.txt.bz2" % (config.data_directory), "r")
header = depth_file.readline()
header = header.strip().split('\t')

#print(len(header))

for line in depth_file:
    line = line.strip().split('\t')
    species_name = line[0]

    if species_name in species_list:
        samples_species = samples_dict[species_name]
        abundance_samples_zip = zip(line[1:], header[1:])
        abundance = [float(x[0]) for x in abundance_samples_zip if x[1] in samples_species]
        mean_abundance = numpy.mean(abundance)

        print(species_name, mean_abundance)

    #print(line)

depth_file.close()
