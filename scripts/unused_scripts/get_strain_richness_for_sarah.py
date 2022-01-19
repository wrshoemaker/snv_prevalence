from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path
import os

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
import scipy.stats as stats
import scipy.special as special

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

#import calculate_predicted_prevalence_mapgd
#import calculate_predicted_prevalence

#import plot_utils


species_list = ['Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', \
        'Ruminococcus_bromii_62047', 'Faecalibacterium_cf_62236', 'Ruminococcus_praunitzii_57453', \
        'Ruminococcus_praunitzii_61481', 'Ruminococcus_praunitzii_62201', 'Prevotella_copri_61740', \
        'Roseburia_inulinivorans_61943', 'Ruminococcus_torques_62045']


intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

count_dict = {}

for species_name in species_list:
    count_dict[species_name] = {}
    count_dict[species_name]['all'] = 0
    count_dict[species_name]['strains'] = 0



for filename in os.listdir('%sstrain_data/' % config.data_directory):
    if filename.endswith(".pkl"):

        filepath = '%sstrain_data/%s' % (config.data_directory, filename)

        with open(filepath, 'rb') as handle:
            b = pickle.load(handle)

        species_in_sample = list(b.keys())

        #print(species_in_sample)

        for species_name in species_list:

            if species_name in species_in_sample:

                count_dict[species_name]['all'] += 1

                if len(b[species_name]) > 1:
                    count_dict[species_name]['strains'] += 1


for species_name in species_list:

    if count_dict[species_name]['all'] > 0:

        fraction_strains = count_dict[species_name]['strains']/count_dict[species_name]['all']

        print(species_name, fraction_strains)
