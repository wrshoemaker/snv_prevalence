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

import matplotlib.pyplot as plt
from scipy.stats import gamma
import scipy.special
from scipy.integrate import quad



good_species_list = parse_midas_data.parse_good_species_list()



for species_name in good_species_list:


    snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

    if os.path.isfile(snp_file_path) == False:
        continue

    # Open post-processed MIDAS output
    snp_file =  bz2.BZ2File(snp_file_path, "r")
    line = snp_file.readline() # header
    items = line.split()[1:]

    samples = numpy.array([item.strip() for item in items])

    bam_file_path = '/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/%s/%_sorted.bam'


Alistipes_finegoldii_56071








'Alistipes_finegoldii_56071' 'Alistipes_onderdonkii_55464' 'Alistipes_putredinis_61533' 'Alistipes_shahii_62199' 'Bacteroidales_bacterium_58650' 'Bacteroides_caccae_53434' 'Bacteroides_cellulosilyticus_58046' 'Bacteroides_fragilis_54507' 'Bacteroides_ovatus_58035' 'Bacteroides_stercoris_56735' 'Bacteroides_thetaiotaomicron_56941' 'Bacteroides_uniformis_57318' 'Bacteroides_vulgatus_57955' 'Bacteroides_xylanisolvens_57185' 'Barnesiella_intestinihominis_62208' 'Dialister_invisus_61905' 'Eubacterium_rectale_56927' 'Oscillibacter_sp_60799' 'Parabacteroides_distasonis_56985' 'Parabacteroides_merdae_56972' 'Ruminococcus_bicirculans_59300' 'Ruminococcus_bromii_62047')
