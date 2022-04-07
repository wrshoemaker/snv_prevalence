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


def get_samples_in_analysis(species_name, clade_type, variant_type):

    frequency_dict = calculate_predicted_prevalence_mapgd.load_frequency_dict(species_name, clade_type)
    samples = []
    for key, value in frequency_dict[variant_type].items():
        samples.extend(value['samples'])

    return list(set(samples))



def get_strain_abundances(species_name, samples):

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    samples_to_keep = []
    richness_to_keep = []
    for sample in samples:

        intermediate_strain_filename = intermediate_strain_filename_template % sample

        if os.path.isfile(intermediate_strain_filename) == False:
            continue

        with open(intermediate_strain_filename, 'rb') as handle:
            b = pickle.load(handle)

        if species_name in b:

            abundances = b[species_name]
            abundances = numpy.asarray(abundances)

            samples_to_keep.append(sample)
            richness_to_keep.append(len(abundances))

    return samples_to_keep, richness_to_keep


#species_name = 'Bacteroides_vulgatus_57955'

#samples = get_samples_in_analysis(species_name, 'all', '4D')

#samples = ['700105882', '700107566c', '700023634c', '700021876', '700038167', '700014498', '700161429', '700014386', '700097482', '700015922c', '700111439', '700034926', '700023390', '700163296', '700105613', '700023845', '700037738', '700095213', '700116865', '700117537', '700037303', '700023788', '700013672', '700037008', '700097234', '700015415', '700015981', '700037967', '700096047', '700033502c', '700016210', '700021361', '700021261', '700117239', '700113867', '700117031', '700116730', #'700023509', '700038806c', '700037649', '700122408', '700037090', '700161967', '700095647', '700038761', '700105049', '700038386', '700024141', '700111156', '700038414', '700106946', '700033363', '700100022', '700164354', '700116917', '700038493', '700164450', '700037405', '700116505', '700021123', '700033435', '700033436', '700015857', '700024437', '700023337', '700024789', '700124561', '700024234', '700024233', '700105991', '700116028', '700161872', '700021295', '700163419', '700021824', #'700034081', '700116296', '700113705', '700095831', '700038956', '700103558', '700109921', '700172813', '700105787', '700111223c', '700099512', '700106170', '700175331', '700172621', '700024711', '700098669', '700106465', '700016765', '700107189', '700034565', '700097688', '700097906c', '700037757', '700165224', '700037357', '700038284', '700161219', '700038870c', '700037852', '700123135', '700099886', '700035785', '700033665', '700161743', '700013603', '700105940', '700037916', '700113810', #'700035373', '700034166', '700032068', '700117069', '700033797', '700021306c', '700097973', '700016610c', '700124429', '700034794', '700024318', '700016000', '700016902', '700023603', '700161342', '700013715', '700024024c', '700109987', '700016960', '700035533', '700038053', '700037284', '700106056', '700103479', '700164249', '700121639', '700167759', '700117766', '700023113', '700014954', '700102043', '700032944', '700105153', '700034254', '700117172', '700105312', '700166862', '700023267', #'700095717', '700105000', '700033201', '700024673', '700106198', '700015113c', '700024349', '700122860', '700098429', '700023578', '700109449', '700033922c', '700034838', '700112376', '700171679', '700024866', '700173546', '700114911', '700117395', '700106754', '700106876c', '700035157', '700023066', '700116668', '700102299', '700023188', '700035861', '700113177', '700033153c', '700023919c', '700013588', '700122765', '700095058', '700013693', '700098867', '700111505', '700033989', '700172069', #'700096700', '700161533', '700098932', '700117303', '700014562c', '700014486', '700016142c', '700102659', '700016542c', '700111745', '700113762', '700023872c', '700172279', '700024752', '700013597', '700024615', '700016456', '700024998', '700106333', '700114218', '700165407', '700175181', '700096380', '700116148']



good_species_list = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                    'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                    'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035',
                    'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318',
                    'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208',
                    'Dialister_invisus_61905', 'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985',
                    'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']


with open(config.data_directory+"strain_dict_sarah.pickle", 'rb') as handle:
    b = pickle.load(handle)


for species_name in good_species_list:

    if species_name not in b:
        print(species_name)
        continue


    b_species = b[species_name]['North America']

    samples = get_samples_in_analysis(species_name, 'all', '4D')

    samples_to_keep = []
    richnes_strains_sarah = []
    for sample in samples:
        if sample in b_species:
            samples_to_keep.append(sample)
            richnes_strains_sarah.append(len(b_species[sample]))



    samples_strain, richnes_strains = get_strain_abundances(species_name, samples_to_keep)

    richnes_strains_sarah = numpy.asarray(richnes_strains_sarah)
    richnes_strains = numpy.asarray(richnes_strains)


    if (len(richnes_strains) == 0) or (len(richnes_strains_sarah) == 0):
        continue

    #print(species_name, sum(richnes_strains>1)/len(richnes_strains), sum(richnes_strains_sarah>1)/len(richnes_strains_sarah))



#print(richnes_strains)
