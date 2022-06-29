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



intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])
#clade_types = ['all','largest_clade', 'no_strains']
clade_types = ['all','largest_clade', 'no_strains', 'just_strains']
#clade_types = ['just_strains']

low_divergence_threshold = config.between_low_divergence_threshold

D_min=20
min_n_samples = 20
min_sample_size = config.between_host_min_sample_size
#min_f_gamma = 0.1
#max_f_gamma = 0.9


#good_species_list = parse_midas_data.parse_good_species_list()
good_species_list = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                    'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                    'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035',
                    'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318',
                    'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208',
                    'Dialister_invisus_61905', 'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985',
                    'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']




f_mean_min_dict = {'Alistipes_finegoldii_56071':0.0005914109589041099, 'Alistipes_onderdonkii_55464':0.00018111538461538463,
                    'Alistipes_putredinis_61533':6.592490118577074e-05, 'Alistipes_shahii_62199':0.0004083333333333334,
                    'Bacteroidales_bacterium_58650':0.0016245517241379325, 'Bacteroides_caccae_53434':0.00016344915254237335,
                    'Bacteroides_cellulosilyticus_58046':0.0005858059701492534, 'Bacteroides_fragilis_54507':0.0005237450980392164,
                    'Bacteroides_ovatus_58035': 5.742696629213457e-05, 'Bacteroides_stercoris_56735': 0.0001016433566433564,
                    'Bacteroides_thetaiotaomicron_56941': 0.00022198333333333376, 'Bacteroides_uniformis_57318':7.263745019920329e-05,
                    'Bacteroides_vulgatus_57955':2.304437869822499e-05, 'Bacteroides_xylanisolvens_57185':0.00012021348314606757,
                    'Barnesiella_intestinihominis_62208':0.0011124038461538467, 'Dialister_invisus_61905': 0.00043426530612244996,
                    'Eubacterium_rectale_56927': 0.0002498376068376069, 'Oscillibacter_sp_60799':0.000893618556701031,
                    'Parabacteroides_distasonis_56985':0.00032050980392156914, 'Parabacteroides_merdae_56972':0.0002769557522123893,
                    'Ruminococcus_bicirculans_59300':9.011688311688348e-05, 'Ruminococcus_bromii_62047':0.0008989999999999998}



# get samples.
# parse mapgd output for all samples
# calculate pi for each sample
# polarize each site using using midas output
# get f_max
# use coverage information from midas
# calculate predicted vs. observed occupancy




def get_samples_by_strain_status(species_name, samples, strain_status=True):

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    s_to_keep = []

    for sample in samples:

        intermediate_strain_filename = intermediate_strain_filename_template % sample

        if os.path.isfile(intermediate_strain_filename) == False:
            continue

        with open(intermediate_strain_filename, 'rb') as handle:
            b = pickle.load(handle)

        if species_name in b:

            abundances = b[species_name]
            # no strains
            if strain_status == False:
                if len(abundances) == 1:
                    s_to_keep.append(str(sample))

            else:
                if len(abundances) > 1:
                    s_to_keep.append(str(sample))


    s_to_keep = numpy.asarray(s_to_keep)

    return s_to_keep



def get_gene_location_name_dict(species_name):

    genome_features_file_path = '%smidas_db_v1.2/rep_genomes/%s/genome.features.gz' % (config.data_directory, species_name)


    gene_chromosome_location_name_dict = {}

    genome_features_file = gzip.open(genome_features_file_path, "r")
    line = genome_features_file.readline() # header
    for line in genome_features_file:
        items = line.strip().split()
        gene_name = items[0]
        chromosome = items[1]
        start = int(items[2])
        stop = int(items[3])
        strand = items[4]
        gene_type = items[5]
        #functions = items[6]

        if chromosome not in gene_chromosome_location_name_dict:
            gene_chromosome_location_name_dict[chromosome] = {}

        gene_chromosome_location_name_dict[chromosome][start] = {}
        gene_chromosome_location_name_dict[chromosome][start]['stop'] = stop
        gene_chromosome_location_name_dict[chromosome][start]['gene_name'] = gene_name

    genome_features_file.close()

    return gene_chromosome_location_name_dict


def get_mapgd_data(species_name, samples, whitelisted_genes):

    gene_location_name_dict = get_gene_location_name_dict(species_name)

    mapgd_dict = {}

    # parse MAPGD output
    mapgd_samples = []
    chromosome_location_tuples = []
    for sample in samples:

        #mapgd_file_path = '/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/%s/%s_sorted.pol' % (species_name, sample)
        #mapgd_file_path = '%smidas_output_bam/%s/%s_sorted.pol' % (config.data_directory, species_name, sample)
        mapgd_file_path = '%smapgd_pol_files/%s/%s_sorted.pol' % (config.data_directory, species_name, sample)

        if os.path.isfile(mapgd_file_path) == False:
            continue

        mapgd_samples.append(sample)

        mapgd_file = open(mapgd_file_path, "r")
        mapgd_header = mapgd_file.readline() # header
        mapgd_header2 = mapgd_file.readline() # header

        #mapgd_dict[sample] = {}

        for line in mapgd_file:
            items = line.strip().split('\t')
            if items[0] == '@END_TABLE':
                continue
            chromosome = items[0].strip()
            location = int(items[1])
            major = items[3]
            minor = items[4]
            coverage = items[5]
            frequency = float(items[7].split('/')[0])

            if chromosome not in gene_location_name_dict:
                continue

            start_locations = list(gene_location_name_dict[chromosome].keys())
            gene_start = max([v for v in start_locations if v <= location] or [None])

            if gene_start == None:
                continue

            # no reading frame
            if location > gene_location_name_dict[chromosome][gene_start]['stop']:
                continue

            gene_name = gene_location_name_dict[chromosome][gene_start]['gene_name']

            if gene_name not in whitelisted_genes:
                continue

            if gene_name not in mapgd_dict:
                mapgd_dict[gene_name] = {}

            if location not in mapgd_dict[gene_name]:
                mapgd_dict[gene_name][location] = {}
                mapgd_dict[gene_name][location]['samples'] = []
                mapgd_dict[gene_name][location]['major_alleles'] = []
                mapgd_dict[gene_name][location]['minor_alleles'] = []
                mapgd_dict[gene_name][location]['frequencies'] = []
                mapgd_dict[gene_name][location]['coverages'] = []

            mapgd_dict[gene_name][location]['samples'].append(sample)
            mapgd_dict[gene_name][location]['major_alleles'].append(major)
            mapgd_dict[gene_name][location]['minor_alleles'].append(minor)
            mapgd_dict[gene_name][location]['frequencies'].append(frequency)
            mapgd_dict[gene_name][location]['coverages'].append(coverage)

            #if (chromosome, location) not in chromosome_location_tuples:
            chromosome_location_tuples.append((chromosome, location))

        mapgd_file.close()

    chromosome_location_tuples = list(set(chromosome_location_tuples))

    return mapgd_samples, mapgd_dict, chromosome_location_tuples



def get_allele_dict(species_name, mapgd_samples, chromosome_location_tuples):

    snps_alt_allele_file_path = '%ssnps/%s/snps_alt_allele.txt.bz2' % (config.data_directory, species_name)
    snps_alt_allele_file = bz2.BZ2File(snps_alt_allele_file_path, "r")
    line = snps_alt_allele_file.readline() # header
    items = line.split()[1:]

    samples = numpy.array([item.strip() for item in items])
    mapgd_samples_idx = numpy.asarray([numpy.where(samples == sample)[0][0] for sample in mapgd_samples])
    allele_dict = {}
    for line in snps_alt_allele_file:
        items = line.split()
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = int(info_items[1])
        ref_allele = info_items[2]

        if (chromosome, location) not in chromosome_location_tuples:
            continue

        alt_alleles = [item for item in items[1:]]
        alt_alleles = numpy.array(alt_alleles)
        alt_alleles = alt_alleles[mapgd_samples_idx]
        alt_alleles = numpy.where(alt_alleles == 'NA', ref_allele, alt_alleles)

        allele_dict[(chromosome, location)] = {}
        allele_dict[(chromosome, location)]['alt_alleles'] = alt_alleles.tolist()
        allele_dict[(chromosome, location)]['ref_allele'] = ref_allele

    snps_alt_allele_file.close()

    return allele_dict



def make_mapgd_alleles_dict(species_to_run='all'):

    if species_to_run == 'all':
        species_to_run = good_species_list

    else:
        species_to_run = [species_to_run]

    intermediate_filename_template = config.data_directory+"mapgd_alleles/%s.dat"

    subject_sample_map = parse_HMP_data.parse_subject_sample_map()

    for species_name in species_to_run:

        sys.stderr.write("%s\n" % species_name)

        gene_location_name_dict = get_gene_location_name_dict(species_name)

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
        snp_file.close()

        samples = numpy.array([item.strip() for item in items])

        # try making the map with this
        samples_all = samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=samples.tolist())]

        sys.stderr.write("Loading MAPGD data...\n")
        #mapgd_samples, mapgd_dict, chromosome_location_tuples = get_mapgd_data(species_name, samples, core_genes)
        mapgd_samples, mapgd_dict, chromosome_location_tuples = get_mapgd_data(species_name, samples_all, core_genes)
        mapgd_samples = numpy.asarray(mapgd_samples)
        #if len(mapgd_samples) == 0:
        #    continue

        #mapgd_samples_idx = numpy.asarray([numpy.where(samples == sample)[0][0] for sample in mapgd_samples])
        # get allelic states
        sys.stderr.write("Loading allelic states...\n")
        allele_dict = get_allele_dict(species_name, mapgd_samples, chromosome_location_tuples)

        intermediate_filename = intermediate_filename_template % species_name

        sys.stderr.write("Saving allelic dict...\n")
        with open(intermediate_filename, 'wb') as handle:
            pickle.dump(allele_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_mapgd_alleles_dict(species_name):

    intermediate_filename_template = config.data_directory+"mapgd_alleles/%s.dat"
    intermediate_filename = intermediate_filename_template % species_name

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)

    return b



def load_pi_dict(species_name):

    intermediate_filename_template = config.data_directory+"pi_annotated_snps_mapgd/%s.dat"
    intermediate_filename = intermediate_filename_template % species_name

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b


def load_parameter_dict(species_name, clade_type):

    intermediate_filename_template = config.data_directory+"parameter_dict_mapgd/%s_%s.dat"
    intermediate_filename = intermediate_filename_template % (species_name, clade_type)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b




def load_abc_status_dict(species_name):

    intermediate_filename_template = config.data_directory+"abc_svm_status_dict/%s.dat"
    intermediate_filename = intermediate_filename_template % species_name

    with open(intermediate_filename, 'rb') as handle:
        #b = pickle.load(handle, protocol=2)
        b = pickle.load(handle)
    return b



def predict_prevalence(f_max, pi_filter, depths, f, min_minor_allele_cov=10):
    # convert to float
    depths = depths.astype(float)
    f_max = float(f_max)

    #pi_filter_non_zero_pi = pi_filter[pi_filter>float(0)]
    #depths_non_zero_pi = depths[pi_filter>float(0)]
    #f_non_zero_pi = f[pi_filter>float(0)]

    pi_filter_non_zero_pi = pi_filter
    depths_non_zero_pi = depths
    f_non_zero_pi = f

    prob_zero = (1 +  f_max * depths_non_zero_pi) ** (-1*pi_filter_non_zero_pi)

    if min_minor_allele_cov > 1:
        # range() excludes upper bound
        cumulative_sum = 0
        for n_i in range(0, min_minor_allele_cov):
            #cumulative_sum += (scipy.special.gamma(pi_filter_non_zero_pi+n_i) / numpy.math.factorial(n_i)) * ((depths_non_zero_pi / (depths_non_zero_pi + (f_max**-1))) ** n_i) / scipy.special.gamma(pi_filter_non_zero_pi)
            cumulative_sum += (scipy.special.gamma(pi_filter_non_zero_pi+n_i) / numpy.math.factorial(n_i)) * ((depths_non_zero_pi / (depths_non_zero_pi + (1/f_max))) ** n_i) / scipy.special.gamma(pi_filter_non_zero_pi)

            #print((depths_non_zero_pi / (depths_non_zero_pi + (f_max**-1)) ) )
            #print( (depths_non_zero_pi / (depths_non_zero_pi + (f_max**-1)) ) ** n_i )

        prob_zero *= cumulative_sum

    predicted_prevalence = 1 - numpy.mean(prob_zero)
    observed_prevalence = sum(f_non_zero_pi>0)/len(f_non_zero_pi)

    return predicted_prevalence, observed_prevalence




def predict_prevalence_slm(freqs, depths, f_mean=None, beta=None, min_minor_allele_cov=10):
    # min_cov=1,
    # min_minor_allele_cov = the coverage where an allele is considered PRESENT
    # for cumulative probability, integral  0 <= n <= (min_minor_allele_cov-1)
    # , folded=False

    # zeros are uninformative here.....
    #freqs = freqs[depths>=min_cov]
    #depths = depths[depths>=min_cov]

    if f_mean == None:
        f_mean = numpy.mean(freqs)

    if beta == None:
        beta = (numpy.mean(freqs)**2)/numpy.var(freqs)


    #if folded == True:
        # prevalence == SEGREGATING
    #    observed_prevalence = sum((freqs>0)&(freqs<1))/len(freqs)

    #else:
    observed_prevalence = sum(freqs>0)/len(freqs)

    prob_zero = (1+ ((f_mean/beta)*depths))**(-1*beta )

    if min_minor_allele_cov > 1:
        # range() excludes upper bound
        cumulative_sum = 0
        # min_minor_allele_cov = upper bound, excluded
        for n_i in range(0, min_minor_allele_cov):
            cumulative_sum += (scipy.special.gamma(beta+n_i) / numpy.math.factorial(n_i)) * ((f_mean*depths/ (beta + f_mean*depths)) ** n_i) / scipy.special.gamma(beta)

        prob_zero *= cumulative_sum


    predicted_prevalence_slm = 1 - numpy.mean(prob_zero)


    return predicted_prevalence_slm, observed_prevalence, f_mean, beta



#def load_predicted_prevalence_subsample_dict():

#    intermediate_filename = config.data_directory+"predicted_prevalence_mapgd_test.dat"

#    with open(intermediate_filename, 'rb') as handle:
#        b = pickle.load(handle)
#    return b



def joint_gamma_parameter_simulation(x_mean_observed, beta_observed, D, distance='euclidean', min_minor_allele_cov=10, iter = 1000, species_name=None):

    if species_name != None:
        x_mean_log10_min = numpy.log10(f_mean_min_dict[species_name])
    else:
        x_mean_log10_min = -4

    n_hosts = len(D)

    x_mean_observed_log10 = numpy.log10(x_mean_observed)
    beta_observed_log10 = numpy.log10(beta_observed)

    x_mean_all_log10 = numpy.random.uniform(low=x_mean_log10_min, high=-0.3, size=iter*100)
    beta_all_log10 = numpy.random.uniform(low=-3, high=2, size=iter*100)

    x_mean_all = 10**x_mean_all_log10
    beta_all = 10**beta_all_log10

    x_mean_truncated_all = []
    beta_truncated_all = []

    x_mean_all_to_keep = []
    beta_all_to_keep = []

    x_mean_prior = []
    beta_prior = []

    n_successes = 0
    n_trials = 0
    #for i in range(iter):
    while n_successes < iter:

        #x_mean_i = x_mean_all[i]
        #beta_i = beta_all[i]
        # get x_mean and betas
        x_mean_i = x_mean_all[n_trials]
        beta_i = beta_all[n_trials]
        # add one to number of trials
        n_trials += 1
        # keep sample sizes the same by generating many samples, truncating them, and selecting n_hosts of them
        #x_i_too_many = stats.gamma.rvs(beta_i, scale=x_mean_i/beta_i, size=n_hosts*1000)
        # numpy is faster
        x_i_too_many = numpy.random.gamma(beta_i, scale=x_mean_i/beta_i, size=n_hosts*10)
        # truncate, empirical allele frequencies are all < 1
        x_i_too_many = x_i_too_many[x_i_too_many<1]

        if len(x_i_too_many) < n_hosts:
            continue

        x_i_to_keep = x_i_too_many[:n_hosts]
        A_i = numpy.random.binomial(D, x_i_to_keep)
        x_prior_i = A_i/D
        x_mean_prior_i = numpy.mean(x_prior_i)
        x_var_prior_i = numpy.var(x_prior_i)
        # we only examine empirical estimates of the variance that are non-zero

        if (x_mean_prior_i == 0) or (x_var_prior_i == 0):
            continue

        beta_prior_i = (x_mean_prior_i**2) / numpy.var(x_prior_i)
        # truncate read counts
        A_i[A_i<min_minor_allele_cov] = 0
        x_truncated_i = A_i/D
        x_mean_truncated_i = numpy.mean(x_truncated_i)
        x_var_truncated_i = numpy.var(x_truncated_i)

        if (x_mean_truncated_i == 0) or (x_var_truncated_i == 0):
            continue

        beta_truncated_i = (x_mean_truncated_i**2) / x_var_truncated_i
        # keep estimates from the prior
        x_mean_prior.append(x_mean_prior_i)
        beta_prior.append(beta_prior_i)

        x_mean_truncated_all.append(x_mean_truncated_i)
        beta_truncated_all.append(beta_truncated_i)

        x_mean_all_to_keep.append(x_mean_i)
        beta_all_to_keep.append(beta_i)

        n_successes += 1

    # was originally "numpy.absolute()". why?
    x_mean_truncated_all = numpy.asarray(x_mean_truncated_all)
    beta_truncated_all = numpy.asarray(beta_truncated_all)

    x_mean_truncated_all_log10 = numpy.log10(x_mean_truncated_all)
    beta_truncated_all_log10 = numpy.log10(beta_truncated_all)

    # get standard deviation of log because you're computing distances on log transformed estimates
    std_dev_x_mean_truncated_log10 = numpy.std(x_mean_truncated_all_log10)
    std_dev_beta_truncated_log10 = numpy.std(beta_truncated_all_log10)

    # euclidean
    #distance_log10 = numpy.sqrt(((x_mean_truncated_all_log10 - x_mean_observed_log10)**2) + ((beta_truncated_all_log10 - beta_observed_log10)**2))
    # normalize by standard deviation

    if distance == 'euclidean':

        distance_log10 = numpy.sqrt((((x_mean_truncated_all_log10 - x_mean_observed_log10)/std_dev_x_mean_truncated_log10)**2) + (((beta_truncated_all_log10 - beta_observed_log10)/std_dev_beta_truncated_log10)**2))

    elif distance == 'manhattan':

        distance_log10 = (numpy.absolute(x_mean_truncated_all_log10 - x_mean_observed_log10)/std_dev_x_mean_truncated_log10) + (numpy.absolute(beta_truncated_all_log10 - beta_observed_log10)/std_dev_beta_truncated_log10)

    # try mean relative error
    #distance_log10 = ((numpy.absolute(x_mean_truncated_all_log10 - x_mean_observed_log10)/x_mean_observed_log10) + (numpy.absolute(beta_truncated_all_log10 - beta_observed_log10)/beta_observed_log10)) / 2

    min_distance = numpy.amin(distance_log10)

    idx_ = numpy.where(distance_log10 == min_distance)[0][0]

    #real_x_mean = x_mean_all_to_keep[idx_]
    #real_beta = beta_all_to_keep[idx_]

    real_x_mean = x_mean_prior[idx_]
    real_beta = beta_prior[idx_]

    return real_x_mean, real_beta, min_distance







def make_pi_and_parameter_dict_and_predict_prevalence(species_to_run='all'):

    subject_sample_map = parse_HMP_data.parse_subject_sample_map()

    if species_to_run == 'all':

        species_to_run = good_species_list

    else:

        species_to_run = [species_to_run]

    #for species_name in good_species_list:
    for species_name in species_to_run:

        sys.stderr.write("%s\n" % species_name)

        gene_location_name_dict = get_gene_location_name_dict(species_name)

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
        # trying this...
        samples_all_unique = samples_all[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=samples_all)]

        sys.stderr.write("Loading MAPGD data...\n")
        #mapgd_samples, mapgd_dict, chromosome_location_tuples = get_mapgd_data(species_name, samples_all, core_genes)
        mapgd_samples, mapgd_dict, chromosome_location_tuples = get_mapgd_data(species_name, samples_all_unique, core_genes)
        # unique samples
        mapgd_samples = numpy.asarray(mapgd_samples)
        mapgd_samples_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples])

        # get samples from largest clade
        snp_samples = diversity_utils.calculate_haploid_samples(species_name)
        if len(snp_samples) < min_sample_size:
            continue
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        # get clade idxs
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = numpy.array(dummy_samples)
        substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        coarse_grained_idxs, coarse_grained_cluster_list = clade_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold)
        coarse_grained_samples = snp_samples[coarse_grained_idxs]
        clade_sets = clade_utils.load_manual_clades(species_name)
        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)
        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])
        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]
        largest_clade_mapgd_samples = list(set(largest_clade_samples) & set(mapgd_samples))
        largest_clade_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in largest_clade_mapgd_samples])

        # no strains
        #mapgd_samples_no_strains = get_samples_no_strains(species_name, mapgd_samples)
        #mapgd_samples_no_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_no_strains])
        mapgd_samples_no_strains = get_samples_by_strain_status(species_name, mapgd_samples, strain_status=False)
        mapgd_samples_no_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_no_strains])


        # just strains
        #mapgd_samples_just_strains = get_samples_just_strains(species_name, mapgd_samples)
        #mapgd_samples_just_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_just_strains])
        mapgd_samples_just_strains = get_samples_by_strain_status(species_name, mapgd_samples, strain_status=True)
        mapgd_samples_just_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_just_strains])


        # get allelic states
        sys.stderr.write("Loading allelic states...\n")
        #allele_dict = get_allele_dict(species_name, mapgd_samples, chromosome_location_tuples)
        allele_dict = load_mapgd_alleles_dict(species_name)

        # get ABC status
        #sys.stderr.write("Loading ABC status dict...\n")
        #abc_status_dict = load_abc_status_dict(species_name)

        sys.stderr.write("Calculating SNP prevalences...\n")
        num_sites_processed = 0

        pi_dict = {}
        for allowed_variant_type in allowed_variant_types:
            pi_dict[allowed_variant_type] = {}
            for sample in mapgd_samples:
                pi_dict[allowed_variant_type][sample] = {}

                pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] = 0
                pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] = 0

                pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] = 0
                pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] = 0

        parameter_dict = {}
        for clade_type in clade_types:
            parameter_dict[clade_type] = {}

        # get chromosome and locations so you can go back and count alt alleles....
        for line in snp_file:

            num_sites_processed+=1
            if num_sites_processed%50000==0:
                sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))
                #if debug:
                #    break

            #if line > 300000:
            #    continue


            items = line.split()
            # Load information about site
            info_items = items[0].split("|")
            chromosome = info_items[0]
            location = long(info_items[1])
            gene_name = info_items[2]
            location_aa = info_items[3]

            # only examine core genes
            if gene_name not in core_genes:
                continue

            if location_aa not in allowed_variant_types:
                continue

            if gene_name not in mapgd_dict:
                continue

            # Load alt and depth counts
            alts = []
            depths = []
            for item in items[1:]:
                subitems = item.split(",")
                alts.append(long(subitems[0]))
                depths.append(long(subitems[1]))
            alts = numpy.array(alts)
            depths = numpy.array(depths)
            refs = depths-alts

            # only try and predict prevalence if there's non-zero coverage for >= 20 hosts
            if sum(depths>=D_min) < min_n_samples:
                continue

            # samples where mapgd was run and

            # only look at sites with mapgd data
            if location not in mapgd_dict[gene_name]:
                continue

            ref_allele = allele_dict[(chromosome, location)]['ref_allele']

            samples = mapgd_dict[gene_name][location]['samples']
            major_alleles = mapgd_dict[gene_name][location]['major_alleles']
            minor_alleles = mapgd_dict[gene_name][location]['minor_alleles']
            # frequencies are of MAJOR alleles
            frequencies = mapgd_dict[gene_name][location]['frequencies']
            coverages = mapgd_dict[gene_name][location]['coverages']

            #major_allele_counts = sorted(major_alleles, key=major_alleles.get, reverse=True)
            # choose the most common major allele as the reference
            frequencies_alts = []
            frequencies_alts_samples = []
            n_samples_alternate_major = 0
            for major_allele_idx, major_allele in enumerate(major_alleles):

                frequency = frequencies[major_allele_idx]
                coverage = coverages[major_allele_idx]
                # we want the frequencies of non reference alleles
                if major_allele == ref_allele:
                    frequency = 1 - frequency
                # count alternative alleles that are also major alles
                else:
                    n_samples_alternate_major += 1

                frequencies_alts.append(frequency)
                frequencies_alts_samples.append(samples[major_allele_idx])

                # only estimate pi at sites with sufficient coverage
                if coverage >= D_min:

                    # add to pi, since 0 < frequency < 1
                    pi_dict[location_aa][samples[major_allele_idx]]['pi_sum_include_boundary'] += 2*frequency*(1-frequency)
                    pi_dict[location_aa][samples[major_allele_idx]]['pi_sum_exclude_boundary'] += 2*frequency*(1-frequency)
                    pi_dict[location_aa][samples[major_allele_idx]]['n_sites_include_boundary'] += 1
                    pi_dict[location_aa][samples[major_allele_idx]]['n_sites_exclude_boundary'] += 1


            for clade_type in clade_types:

                if clade_type == 'largest_clade':
                    alts_iter = alts[largest_clade_idx]
                    depths_iter = depths[largest_clade_idx]
                    refs_iter = refs[largest_clade_idx]
                    samples_to_save = samples_all[largest_clade_idx]
                    samples_clade_type = largest_clade_mapgd_samples


                elif clade_type == 'no_strains':
                    alts_iter = alts[mapgd_samples_no_strains_idx]
                    depths_iter = depths[mapgd_samples_no_strains_idx]
                    refs_iter = refs[mapgd_samples_no_strains_idx]
                    samples_to_save = samples_all[mapgd_samples_no_strains_idx]
                    samples_clade_type = mapgd_samples_no_strains


                elif clade_type == 'just_strains':
                    alts_iter = alts[mapgd_samples_just_strains_idx]
                    depths_iter = depths[mapgd_samples_just_strains_idx]
                    refs_iter = refs[mapgd_samples_just_strains_idx]
                    samples_to_save = samples_all[mapgd_samples_just_strains_idx]
                    samples_clade_type = mapgd_samples_just_strains


                else:
                    alts_iter = alts[mapgd_samples_idx]
                    depths_iter = depths[mapgd_samples_idx]
                    refs_iter = refs[mapgd_samples_idx]
                    samples_to_save = samples_all[mapgd_samples_idx]
                    samples_clade_type = mapgd_samples

                # get samples where MAPGD was run that are also in the clade type

                # original, overwrote by mistake
                #samples = [s for s in samples if s in samples_clade_type]
                #samples_intersection = [s for s in samples if s in samples_clade_type]

                # remove low coverage samples
                alts_iter = alts_iter[depths_iter>=20]
                refs_iter = refs_iter[depths_iter>=20]
                samples_to_save = samples_to_save[depths_iter>=20]
                depths_iter = depths_iter[depths_iter>=20]
                # too few samples, skip
                if len(depths_iter) < min_n_samples:
                    continue

                # get the mapgd samples that are in clade_type
                samples_intersection = [s for s in samples if s in samples_to_save]
                # get non zero frequencies for MAPGD samples that are in the clade type
                frequencies_alts_clade_type = []
                frequencies_alts_clade_type_samples = []
                for frequency_idx, frequency in enumerate(frequencies_alts):
                    # keep the frequency and get the sample if its in the clade_type
                    if frequencies_alts_samples[frequency_idx] in samples_intersection:
                        frequencies_alts_clade_type.append(frequency)
                        frequencies_alts_clade_type_samples.append(frequencies_alts_samples[frequency_idx])

                # number of zeros
                #samples_with_absence = set(samples_clade_type) - set(samples_intersection)
                #samples_with_absence = set(samples_clade_type) - set(frequencies_alts_clade_type_samples)
                # get mapgd samples that are in your clade_type that didn't have non-zero frequencies
                samples_with_absence = set(samples_to_save) - set(frequencies_alts_clade_type_samples)

                #n_zeros_to_add = len(mapgd_samples) - len(samples)
                # add those zeros back in!
                frequencies_alts_clade_type.extend([0]*len(samples_with_absence))
                frequencies_alts_clade_type = numpy.asarray(frequencies_alts_clade_type)
                # add samples back
                frequencies_alts_clade_type_samples.extend(samples_with_absence)

                # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
                #population_prevalence = n_samples_alternate_major/len(samples_clade_type)
                population_prevalence = n_samples_alternate_major/len(major_alleles)
                reference = True
                if population_prevalence>0.5:
                    frequencies_alts_clade_type = 1 - frequencies_alts_clade_type
                    reference = False

                f_max = max(frequencies_alts_clade_type)

                #if f_max == float(1):
                #    continue

                # ignore sites where all hosts are fixed for the alternate allele
                if sum(frequencies_alts_clade_type == float(1)) == len(frequencies_alts_clade_type):
                    continue

                f_mean = numpy.mean(frequencies_alts_clade_type)
                f_var = numpy.var(frequencies_alts_clade_type)

                prevalence = sum(frequencies_alts_clade_type>0)/len(frequencies_alts_clade_type)

                frequencies_alts_clade_type_no_zeros = frequencies_alts_clade_type[frequencies_alts_clade_type>0]

                if prevalence == 0:
                    continue

                if len(frequencies_alts_clade_type_no_zeros) > 3:
                    f_mean_no_zeros = numpy.mean(frequencies_alts_clade_type_no_zeros)
                    f_var_no_zeros = numpy.var(frequencies_alts_clade_type_no_zeros)
                else:
                    f_mean_no_zeros = float("nan")
                    f_var_no_zeros = float("nan")


                # go through all the samples where the allele was absent
                for s in samples_with_absence:
                    pi_dict[location_aa][s]['n_sites_include_boundary'] += 1

                parameter_dict[clade_type][(chromosome, location)] = {}

                parameter_dict[clade_type][(chromosome, location)]['f_max'] = f_max
                parameter_dict[clade_type][(chromosome, location)]['f_mean'] = f_mean
                parameter_dict[clade_type][(chromosome, location)]['f_var'] = f_var
                parameter_dict[clade_type][(chromosome, location)]['observed_prevalence'] = prevalence
                parameter_dict[clade_type][(chromosome, location)]['depths'] = depths_iter
                #parameter_dict[clade_type][(chromosome, location)]['samples'] = samples_to_save
                parameter_dict[clade_type][(chromosome, location)]['samples'] = frequencies_alts_clade_type_samples
                parameter_dict[clade_type][(chromosome, location)]['location_aa'] = location_aa
                parameter_dict[clade_type][(chromosome, location)]['f_no_zeros_mapgd'] = frequencies_alts_clade_type_no_zeros.tolist()
                parameter_dict[clade_type][(chromosome, location)]['f_mapgd'] = frequencies_alts_clade_type.tolist()
                parameter_dict[clade_type][(chromosome, location)]['f_mean_no_zeros'] = f_mean_no_zeros
                parameter_dict[clade_type][(chromosome, location)]['f_var_no_zeros'] = f_var_no_zeros




        snp_file.close()

        # calculate pi
        for allowed_variant_type in allowed_variant_types:
            for sample in mapgd_samples:
                # only keep a host
                if pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] > 0:
                    pi_dict[allowed_variant_type][sample]['pi_include_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_include_boundary']
                    pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] = pi_dict[allowed_variant_type][sample]['n_sites_include_boundary']

                    #pi_dict[allowed_variant_type][sample]['pi_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']
                    #pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']

                if pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] > 0:
                    pi_dict[allowed_variant_type][sample]['pi_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']
                    pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']


        for clade_type in clade_types:
            prevalence_dict_mapgd = {}
            frequency_dict_mapgd = {}

            #mapgd_samples_to_keep_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples])
            #things_to_measure = ['f_max_all', 'observed_prevalence_all', 'f_mean_all', 'predicted_prevalence_all', 'predicted_f_mean_all', 'f_max_mapgd', 'observed_prevalence_mapgd', 'predicted_prevalence_mapgd', 'f_mean_mapgd', 'predicted_f_mean_mapgd', 'f_no_zeros_mapgd', 'f_var_mapgd', 'predicted_f_var_mapgd', 'n_non_zero_frequencies', 'sites', 'mean_coverage_mapgd', 'predicted_prevalence_mapgd_slm', 'observed_prevalence_mapgd_slm']

            things_to_measure = ['f_max_all', 'observed_prevalence_all', 'f_mean_all', 'predicted_prevalence_all', \
                                'predicted_f_mean_all', 'f_max_mapgd', 'observed_prevalence_mapgd', 'predicted_prevalence_mapgd', \
                                'f_mean_mapgd', 'predicted_f_mean_mapgd', 'f_no_zeros_mapgd', 'f_var_mapgd', 'predicted_f_var_mapgd', \
                                'n_non_zero_frequencies', 'sites', 'mean_coverage_mapgd', 'observed_prevalence_mapgd_slm', 'predicted_prevalence_mapgd_slm',\
                                'f_mean_slm', 'beta_slm', 'f_mean_slm_best', 'beta_slm_best', 'min_distance_slm', 'predicted_prevalence_mapgd_slm_best', 'observed_prevalence_mapgd_slm_best']
                                #'f_max_mapgd_folded', 'predicted_prevalence_mapgd_folded', 'f_mean_mapgd_folded', \
                                #'predicted_f_mean_mapgd_folded', 'f_var_mapgd_folded', 'predicted_f_var_mapgd_folded',\
                                #'f_no_zeros_mapgd_folded', 'predicted_prevalence_mapgd_slm_folded', 'observed_prevalence_mapgd_slm_folded']

            for allowed_variant_type in allowed_variant_types:
                frequency_dict_mapgd[allowed_variant_type] = {}

            #'f_mean_no_zeros_mapgd', 'f_var_no_zeros_mapgd', 'n_hosts'
            # observed_prevalence_mapgd_slm_cov_20,  observed_prevalence_mapgd_slm
            #for pi_type in ['pi_exclude_boundary', 'pi_include_boundary']:
            for pi_type in ['pi_include_boundary']:
                prevalence_dict_mapgd[pi_type] = {}
                for allowed_variant_type in allowed_variant_types:
                    prevalence_dict_mapgd[pi_type][allowed_variant_type] = {}
                    for t in things_to_measure:
                        prevalence_dict_mapgd[pi_type][allowed_variant_type][t] = []

            sys.stderr.write("%s\n" % clade_type)
            sys.stderr.write("Getting parameter dict...\n")
            parameter_dict_clade = parameter_dict[clade_type]
            sys.stderr.write("Done!\n")
            sys.stderr.write("Making prevalence dict...\n")

            for key, value in parameter_dict_clade.items():

                # use the depth of coverage used by MAPGD
                depths_i = value['depths']
                f_i = value['f_mapgd']
                samples_i = value['samples']
                location_aa_i = value['location_aa']
                f_max_i = value['f_max']
                f_mean_i = value['f_mean']
                f_var_i = value['f_var']
                f_no_zeros_i = value['f_no_zeros_mapgd']
                f_mean_no_zeros_i = value['f_mean_no_zeros']
                f_var_no_zeros_i = value['f_var_no_zeros']

                #for pi_type in ['pi_exclude_boundary', 'pi_include_boundary']:
                for pi_type in ['pi_include_boundary']:

                    pi_filter = []
                    depths_i_filter = []
                    f_i_filter = []
                    samples_i_filter = []

                    for sample_idx, sample in enumerate(samples_i):
                        if sample in pi_dict[location_aa_i]:
                            if pi_type in pi_dict[location_aa_i][sample]:
                                #mapgd_samples_to_keep.append(sample)
                                depths_i_idx = depths_i[sample_idx]
                                pi_filter.append(pi_dict[location_aa_i][sample][pi_type])
                                #depths_i_filter.append(depths_i[sample_idx])
                                depths_i_filter.append(depths_i_idx)
                                f_i_filter.append(f_i[sample_idx])
                                samples_i_filter.append(sample)


                    pi_filter = numpy.asarray(pi_filter)
                    depths_i_filter = numpy.asarray(depths_i_filter)
                    f_i_filter = numpy.asarray(f_i_filter)
                    samples_i_filter = numpy.asarray(samples_i_filter)

                    if len(pi_filter) < min_n_samples:
                        continue

                    f_i_filter_non_zero = f_i_filter[f_i_filter>0]
                    samples_i_filter_non_zero = samples_i_filter[f_i_filter>0]
                    #predicted_prevalence = 1 - numpy.mean( (1 +  f_max_i * depths_i_filter) ** (-1*pi_filter) )

                    #print(f_max_i, pi_filter, depths_i_filter)
                    predicted_prevalence, observed_prevalence = predict_prevalence(f_max_i, pi_filter, depths_i_filter, f_i_filter, min_minor_allele_cov=10)

                    predicted_f_mean = f_max_i * numpy.mean(pi_filter)
                    predicted_f_var = (f_max_i**2) * numpy.mean(pi_filter)

                    mean_coverage_mapgd = numpy.mean(depths_i_filter)
                    #observed_prevalence = sum(f_i_filter>0)/len(f_i_filter)
                    predicted_prevalence_slm, observed_prevalence_slm, f_mean_slm, beta_slm = predict_prevalence_slm(f_i_filter, depths_i_filter, min_minor_allele_cov=10)

                    if predicted_prevalence_slm == float('nan'):
                        continue

                    #if clade_type not in abc_status_dict:
                    #    abc_status = False

                    #else:
                    #    if location_aa_i not in abc_status_dict[clade_type]:
                    #        abc_status = False

                    #    else:
                    #        if key in abc_status_dict[clade_type][location_aa_i]:
                    #            abc_status = abc_status_dict[clade_type][location_aa_i][key]
                    #        else:
                    #            abc_status = False


                    #if abc_status == True:
                    #f_mean_best, beta_best, min_distance_slm = joint_gamma_parameter_simulation(f_mean_slm, beta_slm, depths_i_filter, species_name=species_name)
                    #predicted_prevalence_slm_best, observed_prevalence_slm_best, f_mean_slm_best, beta_slm_best = predict_prevalence_slm(f_i_filter, depths_i_filter, f_mean=f_mean_best, beta=beta_best, min_minor_allele_cov=10)

                    f_mean_best = float('nan')
                    beta_best = float('nan')
                    min_distance_slm = float('nan')

                    predicted_prevalence_slm_best = float('nan')
                    observed_prevalence_slm_best = float('nan')
                    f_mean_slm_best = float('nan')
                    beta_slm_best = float('nan')

                    #else:
                    #    f_mean_best = float('nan')
                    #    beta_best = float('nan')
                    #    min_distance_slm = float('nan')

                    #    predicted_prevalence_slm_best = float('nan')
                    #    observed_prevalence_slm_best = float('nan')
                    #    f_mean_slm_best = float('nan')
                    #    beta_slm_best = float('nan')

                    #error = numpy.absolute(observed_prevalence_slm_best - predicted_prevalence_slm_best)/observed_prevalence_slm_best
                    #print(observed_prevalence_slm_best, predicted_prevalence_slm_best, error)

                    # add frequency data.
                    # we dont need to save the zeros
                    if pi_type == 'pi_include_boundary':
                        frequency_dict_mapgd[location_aa_i][key] = {}
                        frequency_dict_mapgd[location_aa_i][key]['frequencies'] = f_i_filter_non_zero.tolist()
                        frequency_dict_mapgd[location_aa_i][key]['samples_non_zero'] = samples_i_filter_non_zero.tolist()
                        frequency_dict_mapgd[location_aa_i][key]['samples'] = samples_i_filter.tolist()
                        frequency_dict_mapgd[location_aa_i][key]['n_samples'] = len(samples_i_filter)


                    prevalence_dict_mapgd[pi_type][location_aa_i]['f_max_mapgd'].append(f_max_i)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['observed_prevalence_mapgd'].append(observed_prevalence)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['predicted_prevalence_mapgd'].append(predicted_prevalence)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['predicted_f_mean_mapgd'].append(predicted_f_mean)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['f_mean_mapgd'].append(f_mean_i)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['f_var_mapgd'].append(f_var_i)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['f_mean_mapgd_no_zeros'].append(f_mean_no_zeros_i)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['f_var_mapgd_no_zeros'].append(f_var_no_zeros_i)

                    prevalence_dict_mapgd[pi_type][location_aa_i]['predicted_f_var_mapgd'].append(predicted_f_var)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['f_no_zeros_mapgd'].extend(f_no_zeros_i)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['n_non_zero_frequencies'].append(len(f_i_filter_non_zero))
                    prevalence_dict_mapgd[pi_type][location_aa_i]['mean_coverage_mapgd'].append(mean_coverage_mapgd)
                    # predicted_prevalence_mapgd_slm
                    prevalence_dict_mapgd[pi_type][location_aa_i]['predicted_prevalence_mapgd_slm'].append(predicted_prevalence_slm)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['observed_prevalence_mapgd_slm'].append(observed_prevalence_slm)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['f_mean_slm'].append(f_mean_slm)
                    prevalence_dict_mapgd[pi_type][location_aa_i]['beta_slm'].append(beta_slm)

                    #prevalence_dict_mapgd[pi_type][location_aa_i]['f_mean_slm_best'].append(f_mean_best)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['beta_slm_best'].append(beta_best)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['min_distance_slm'].append(min_distance_slm)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['predicted_prevalence_mapgd_slm_best'].append(predicted_prevalence_slm_best)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['observed_prevalence_mapgd_slm_best'].append(observed_prevalence_slm_best)

                    #prevalence_dict_mapgd[pi_type][location_aa_i]['predicted_prevalence_slm_minor_cov_threshold'].append(predicted_prevalence_slm_minor_cov_threshold)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['observed_prevalence_slm_minor_cov_threshold'].append(observed_prevalence_slm_minor_cov_threshold)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['f_mean_slm_minor_cov_threshold'].append(f_mean_slm_minor_cov_threshold)
                    #prevalence_dict_mapgd[pi_type][location_aa_i]['f_beta_slm_minor_cov_threshold'].append(f_beta_slm_minor_cov_threshold)

                    prevalence_dict_mapgd[pi_type][location_aa_i]['sites'].append(key)


            sys.stderr.write("Saving prevalence dict...\n")
            intermediate_filename = "%spredicted_prevalence_dict_mapgd/%s_%s.dat" % (config.data_directory, species_name, clade_type)
            with open(intermediate_filename, 'wb') as handle:
                pickle.dump(prevalence_dict_mapgd, handle, protocol=pickle.HIGHEST_PROTOCOL)

            sys.stderr.write("Saving frequency dict...\n")
            intermediate_filename = "%sfrequency_dict_mapgd/%s_%s.dat" % (config.data_directory, species_name, clade_type)
            with open(intermediate_filename, 'wb') as handle:
                pickle.dump(frequency_dict_mapgd, handle, protocol=pickle.HIGHEST_PROTOCOL)



def make_frequency_dict_mapgd_non_zero():

    intermediate_path = "%sfrequency_dict_mapgd/" % (config.data_directory)

    for filename in os.listdir(intermediate_path):
        if filename.endswith(".dat"):
            file_path = os.path.join(intermediate_path, filename)

            sys.stderr.write("%s\n" % filename)

            with open(file_path, 'rb') as handle:
                frequency_dict_mapgd = pickle.load(handle)

            variant_types = list(frequency_dict_mapgd.keys())

            frequency_dict_mapgd_non_zero = {}

            for variant_type in variant_types:

                frequency_dict_mapgd_non_zero[variant_type] = []

                for chromosome_location_tuple in frequency_dict_mapgd[variant_type].keys():

                    freqs_non_zero = [f for f in frequency_dict_mapgd[variant_type][chromosome_location_tuple]['frequencies'] if f > 0]

                    frequency_dict_mapgd_non_zero[variant_type].extend(freqs_non_zero)

            sys.stderr.write("Saving non-zero frequency dict...\n")
            intermediate_filename_non_zero = "%sfrequency_dict_mapgd_non_zero/%s" % (config.data_directory, filename)
            with open(intermediate_filename_non_zero, 'wb') as handle:
                pickle.dump(frequency_dict_mapgd_non_zero, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_frequency_dict_mapgd_non_zero(species_name, clade_type):

    #intermediate_filename = config.data_directory+"predicted_prevalence_mapgd_test.dat"
    intermediate_filename = "%sfrequency_dict_mapgd_non_zero/%s_%s.dat" % (config.data_directory, species_name, clade_type)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b



def load_predicted_prevalence_dict(species_name, clade_type):

    #intermediate_filename = config.data_directory+"predicted_prevalence_mapgd_test.dat"
    intermediate_filename = "%spredicted_prevalence_dict_mapgd/%s_%s.dat" % (config.data_directory, species_name, clade_type)

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b


def load_frequency_dict(species_name, clade_type):

    intermediate_filename = "%sfrequency_dict_mapgd/%s_%s.dat" % (config.data_directory, species_name, clade_type)
    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b



def load_predicted_prevalence_dict_all(test=False):

    intermediate_path = "%spredicted_prevalence_dict_mapgd/" % (config.data_directory)

    prevalence_dict = {}

    for filename in os.listdir(intermediate_path):
        if filename.endswith(".dat"):
            file_path = os.path.join(intermediate_path, filename)
            filename = filename.split('.')[0]

            if filename.count('_') == 3:
                split = filename.rsplit('_', 1)
                species = split[0]
                clade_type = split[1]
            else:
                last_delim = 2
                split = filename.rsplit('_', 2)
                species = split[0]
                clade_type = split[1] + '_' + split[2]

            if test == True:
                #if species != 'Bacteroides_vulgatus_57955':
                if species != 'Bacteroides_vulgatus_57955':
                    continue

            if species not in prevalence_dict:
                prevalence_dict[species] = {}


            with open(file_path, 'rb') as handle:
                b = pickle.load(handle)
            #return b
            prevalence_dict[species][clade_type] = b



    return prevalence_dict







if __name__=='__main__':

    #calculate_pi_from_snps_fmax_cutoff(f_max_range)

    species_to_run = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                        'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                        'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035',
                        'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318',
                        'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208',
                        'Dialister_invisus_61905', 'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985',
                        'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", help="Name of specific species to run code on", default='none')
    #parser.add_argument("--species", help="Name of specific species to run code on", default='none')
    parser.add_argument("--make_alleles_dict", help="Get alleleic states for MAPGD", default='all')
    #parser.add_argument("--make_pi_and_parameter_dict", help="Get pi and parameter dict", default='none')
    #parser.add_argument("--make_frequency_dict_mapgd_non_zero", help="Get non zero frequencies for all sites used in analysis", default='none')
    parser.add_argument("--make_frequency_dict_mapgd_non_zero", help="Get non zero frequencies for all sites used in analysis", action='store_true')



    args = parser.parse_args()

    species=args.species
    make_alleles_dict=args.make_alleles_dict
    #make_pi_and_parameter_dict=args.make_pi_and_parameter_dict
    make_frequency_dict_mapgd_non_zero_ = args.make_frequency_dict_mapgd_non_zero


    if species != 'none':
        #make_pi_and_parameter_dict(species)
        make_pi_and_parameter_dict_and_predict_prevalence(species)


    if make_alleles_dict != 'none':
        make_mapgd_alleles_dict(species_to_run=make_alleles_dict)

    if make_frequency_dict_mapgd_non_zero_ == True:
        make_frequency_dict_mapgd_non_zero()


    #if make_pi_and_parameter_dict != 'none':
    #    make_pi_and_parameter_dict(species)
