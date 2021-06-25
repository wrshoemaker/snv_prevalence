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

import matplotlib.pyplot as plt




intermediate_filename_template = '%s%s.txt.gz'

pairwise_directory = '%spairwise_distance/' % (parse_midas_data.data_directory)


#good_species_list = parse_midas_data.parse_good_species_list()

good_species_list = ['Bacteroides_vulgatus_57955']


for species_name in good_species_list:

    intermediate_filename = intermediate_filename_template % (pairwise_directory, species_name)

    if os.path.isfile(intermediate_filename) == False:
        continue

    file = gzip.open(intermediate_filename,"r")
    header_line = file.readline() # header
    header_items = header_line.split(",")

    sample_pairs = []
    substitution_rates = []
    gene_overlap = []
    for line in file:

        line_split = line.strip().split(', ')

        sample_pairs.append((line_split[0], line_split[1]))
        substitution_rates.append(float(line_split[2]))
        gene_overlap.append(float(line_split[5]))

    substitution_rates = numpy.asarray(substitution_rates)
    gene_overlap = numpy.asarray(gene_overlap)




    fig, ax = plt.subplots(figsize=(4,4))


    ax.scatter(substitution_rates, gene_overlap, alpha=0.3)
    ax.set_xlim([min(substitution_rates), max(substitution_rates)])

    ax.set_xlabel( 'Genetic distance in core genes' )
    ax.set_ylabel('Proportion of shared genes' )

    ax.set_title(species_name)

    #ax.set_xscale("log", basex=10)

    fig.savefig("%sdistance_gene_overlap_%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
