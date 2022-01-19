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
import scipy.stats as stats

import diversity_utils
import core_gene_utils
import parse_midas_data
import parse_HMP_data
import sample_utils
import calculate_substitution_rates
import clade_utils

from scipy.stats import gamma
import scipy.special
from scipy.integrate import quad

import calculate_predicted_prevalence_mapgd


def get_samples_just_strains(species_name, samples):

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    samples_no_strains = []

    for sample in samples:

        intermediate_strain_filename = intermediate_strain_filename_template % sample

        if os.path.isfile(intermediate_strain_filename) == False:
            continue

        with open(intermediate_strain_filename, 'rb') as handle:
            b = pickle.load(handle)

        if species_name in b:

            abundances = b[species_name]

            if len(abundances) > 1:

                samples_no_strains.append(str(sample))

    samples_no_strains = numpy.asarray(samples_no_strains)

    return samples_no_strains



def get_samples_no_strains(species_name, samples):

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    samples_no_strains = []

    for sample in samples:

        intermediate_strain_filename = intermediate_strain_filename_template % sample

        if os.path.isfile(intermediate_strain_filename) == False:
            continue

        with open(intermediate_strain_filename, 'rb') as handle:
            b = pickle.load(handle)

        if species_name in b:

            abundances = b[species_name]

            if len(abundances) == 1:

                samples_no_strains.append(str(sample))

    samples_no_strains = numpy.asarray(samples_no_strains)

    return samples_no_strains





species_list = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533', \
                'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434', \
                'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035', \
                'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318', \
                'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208', \
                'Dialister_invisus_61905',  'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985', \
                'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']



for species_name in species_list:

    sys.stderr.write("%s\n" % species_name)

    gene_location_name_dict = calculate_predicted_prevalence_mapgd.get_gene_location_name_dict(species_name)

    sys.stderr.write("Loading whitelisted genes...\n")
    core_genes = core_gene_utils.parse_core_genes(species_name)

    # Holds panel wide prevalence for each species
    #os.system('mkdir -p %ssnp_prevalences' % config.data_directory)
    snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

    if os.path.isfile(snp_file_path) == False:
        continue

    # Open post-processed MIDAS output
    snp_file =  bz2.BZ2File(snp_file_path, "r")
    line = snp_file.readline() # header
    items = line.split()[1:]

    samples_all = numpy.array([item.strip() for item in items])

    sys.stderr.write("Loading MAPGD data...\n")
    mapgd_samples, mapgd_dict, chromosome_location_tuples = calculate_predicted_prevalence_mapgd.get_mapgd_data(species_name, samples_all, core_genes)


    print(species_name, len(get_samples_just_strains(species_name, mapgd_samples)), len(get_samples_no_strains(species_name, mapgd_samples)))



    snp_file.close()
