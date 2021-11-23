from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

from collections import Counter


import diversity_utils
import figure_utils
import parse_midas_data



import calculate_predicted_prevalence
import calculate_predicted_prevalence_mapgd


prevalence_dict = calculate_predicted_prevalence.load_predicted_prevalence_subsample_dict()
prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_subsample_dict()


species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()

clade_type = 'all'
variant_type = '4D'

proportion_dict = {}
pi_all_dict = {}

for species_name in species_list:

    pi_dict = calculate_predicted_prevalence.load_pi_dict(species_name)
    pi_dict_mapgd = calculate_predicted_prevalence_mapgd.load_pi_dict(species_name)

    proportion_dict[species_name] = {}
    pi_all_dict[species_name] = {}

    n_sampled_mapgd = 0

    for sample, sample_pi_dict in pi_dict[variant_type].items():

        if ('pi_exclude_boundary' in sample_pi_dict) and ('pi_include_boundary' in sample_pi_dict):

            proportion_zeros = 1 - sample_pi_dict['n_sites_exclude_boundary']/sample_pi_dict['n_sites_include_boundary']

            pi_include_boundary = sample_pi_dict['pi_include_boundary']
            pi_exclude_boundary = sample_pi_dict['pi_exclude_boundary']

            if sample not in proportion_dict[species_name]:
                proportion_dict[species_name][sample] = {}

            if sample not in pi_all_dict[species_name]:
                pi_all_dict[species_name][sample] = {}

            pi_all_dict[species_name][sample]['naive'] = {}
            pi_all_dict[species_name][sample]['naive']['pi_include_boundary'] = pi_include_boundary
            pi_all_dict[species_name][sample]['naive']['pi_exclude_boundary'] = pi_exclude_boundary


    for sample, sample_pi_dict in pi_dict_mapgd[variant_type].items():

        if ('pi_exclude_boundary' in sample_pi_dict) and ('pi_include_boundary' in sample_pi_dict):

            proportion_zeros = 1 - sample_pi_dict['n_sites_exclude_boundary']/sample_pi_dict['n_sites_include_boundary']

            if sample not in proportion_dict[species_name]:
                proportion_dict[species_name][sample] = {}

            #proportion_dict[species_name][sample]['mapgd'] = proportion_zeros


            pi_include_boundary = sample_pi_dict['pi_include_boundary']
            pi_exclude_boundary = sample_pi_dict['pi_exclude_boundary']

            if sample not in pi_all_dict[species_name]:
                pi_all_dict[species_name][sample] = {}


            pi_all_dict[species_name][sample]['mapgd'] = {}
            pi_all_dict[species_name][sample]['mapgd']['pi_include_boundary'] = pi_include_boundary
            pi_all_dict[species_name][sample]['mapgd']['pi_exclude_boundary'] = pi_exclude_boundary

            n_sampled_mapgd += 1



    print(figure_utils.get_pretty_species_name(species_name), n_sampled_mapgd)





#print(proportion_dict)
