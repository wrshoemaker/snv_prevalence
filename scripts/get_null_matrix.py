
import sys
import config
import parse_midas_data
import parse_HMP_data
import diversity_utils
import numpy
import sample_utils

from numpy.random import randint, random, choice, multinomial, shuffle



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
parser.add_argument("--modification-threshold", type=int, help="max number of SNV differences before calling a modification", default=config.modification_difference_threshold)



args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize
modification_difference_threshold = args.modification_threshold
replacement_difference_threshold = config.replacement_difference_threshold
twin_modification_difference_threshold = config.twin_modification_difference_threshold
twin_replacement_difference_threshold = config.twin_replacement_difference_threshold
default_num_bootstraps = 10000



################################################################################

#####################
#
# Settings for calculation:
#
#####################

min_coverage = config.min_median_coverage
min_sample_size = 3
min_haploid_sample_size = 10

variant_types = ['1D','4D']

within_host_type = 'consecutive' # consecutive timepoints

# For partitioning SNVs according to prevalence
derived_freq_bins = numpy.array([-1,0,0.01,0.1,0.5,0.9,0.99,1,2])
derived_virtual_freqs = numpy.arange(0,len(derived_freq_bins)-1)
derived_virtual_xticks = list(derived_virtual_freqs[:-1]+0.5)
derived_virtual_xticklabels = ['0','.01','.1','.5','.9','.99','1']

# For partitioning genes into different prevalence classes
gene_freq_bins = numpy.array([-1,0.1,0.5,0.9,2])
gene_freq_xticks      = [-4, -3,  -2,   -1,   0,   1,    2,   3, 4]
gene_freq_xticklabels = ['0','0.1','0.5', '0.9','1','0.9','0.5', '0.1','0']
gene_gain_virtual_freqs = numpy.array([3.5,2.5,1.5,0.5])
gene_loss_virtual_freqs = numpy.array([-3.5,-2.5,-1.5,-0.5])

#####
#
# Settings for different cohorts we are looking at
#
#####
#cohorts = ["hmp", "twins", "young_twins"]
cohorts = ["hmp"]
countries = ["United States", "United Kingdom", "Western Europe"]
country_cohort_map = {country: cohort for country,cohort in zip(countries,cohorts)}

modification_difference_thresholds = {"hmp": modification_difference_threshold, "twins": 1e06, "young_twins": twin_modification_difference_threshold}

replacement_difference_thresholds = {"hmp": replacement_difference_threshold, "twins": 1e06, "young_twins": twin_replacement_difference_threshold}




species_name='Bacteroides_vulgatus_57955'

# Load gene presence/absence information for species_name
#sys.stderr.write("Loading %s...\n" % species_name)
#gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
#sys.stderr.write("Done!\n")
# remove genes that occur in all samples ???



samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name,     allowed_variant_types=set(['4D']))

print(sfs_map[b'700015181'])



################################
#
# Now do calculation
#
################################

hmp_species_qp_counts = {}
twin_species_qp_counts = {}

species_snp_change_distribution = {cohort: {} for cohort in cohorts}
species_snp_nerrs = {cohort: {} for cohort in cohorts}
species_gene_change_distribution = {cohort: {} for cohort in cohorts}
species_gene_nerrs = {cohort: {} for cohort in cohorts}

# observed within host value
pooled_snp_change_distribution = {cohort : [] for cohort in cohorts}
pooled_gene_change_distribution = {cohort : [] for cohort in cohorts} # for modifications

# typical value, median other sample
pooled_between_snp_change_distribution = {cohort : [] for cohort in cohorts}
pooled_between_gene_change_distribution = {cohort : [] for cohort in cohorts}

# closest other sample
pooled_min_between_snp_change_distribution = {cohort : [] for cohort in cohorts}
pooled_min_between_gene_change_distribution = {cohort : [] for cohort in cohorts}

replacement_map = {cohort: {} for cohort in cohorts}

total_freq_snps = {cohort: {} for cohort in cohorts}
total_null_freq_snps = {cohort: {} for cohort in cohorts}
for cohort in cohorts:

    # This holds the # of SNVs in each prevalence class of each var type
    total_freq_snps[cohort] = {var_type: numpy.zeros_like(derived_virtual_freqs) for var_type in variant_types}

    # This holds the expectation of the # of SNVs in each prevalence class of 1D and 4D
    # (i.e., relative fraction of opportunities for different species), conditioned on prevalence
    total_null_freq_snps[cohort] = {var_type: numpy.zeros_like(derived_virtual_freqs) for var_type in variant_types}


# This sums up across the different var_type categories
total_freq_all_snps = {cohort: numpy.zeros_like(derived_virtual_freqs) for cohort in cohorts}

# This is the null distribution of prevalence (conditioned on total # of SNVs)
total_null_freq_all_snps = {cohort: numpy.zeros_like(derived_virtual_freqs)*1.0 for cohort in cohorts}

total_freq_gains = {cohort: numpy.zeros(len(gene_freq_bins)-1)*1.0 for cohort in cohorts}
total_freq_losses = {cohort: numpy.zeros_like(total_freq_gains[cohort]) for cohort in cohorts}
total_null_freq_losses = {cohort: numpy.zeros_like(total_freq_gains[cohort]) for cohort in cohorts}


# SNV and gene prevalences resolved by subject
# so that we can do bootstrapping
snv_prevalence_count_map = {cohort: {} for cohort in cohorts}
gene_gain_count_map = {cohort: {} for cohort in cohorts}
gene_loss_count_map = {cohort: {} for cohort in cohorts}

snv_prevalence_map = {cohort : {} for cohort in cohorts}
gene_gain_prevalence_map = {cohort: {} for cohort in cohorts}
gene_loss_prevalence_map = {cohort: {} for cohort in cohorts}

variant_type_prevalence_map = {cohort: {'1D':[], '4D':[]} for cohort in cohorts}
output_strs = {cohort : [] for cohort in cohorts}





######
#
# Helper function for calculating the prevalence of the sweeping allele in the larger cohort
#
######
def get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map):

    gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change

    f1 = A1*1.0/D2
    f2 = A2*1.0/D2

    is_reversion = (f1>f2)

    location_tuple = (contig, position)

    is_private_snv = (location_tuple in private_snv_map)


    # Now calculate frequency-stratified version

    if location_tuple in snv_freq_map:
        f = snv_freq_map[location_tuple]
    else:
        sys.stderr.write("SNP not in map. Shouldn't happen!\n")
        f = -0.5

    # Let's impose that private snvs have zero freq (specifically, lim 0^-)
    if is_private_snv:
        f = -0.5

    # Change f so that it represents
    # frequency of allele at second timepoint
    if is_reversion:
        f = 1-f

    return f


# Helper functions for stats
from scipy.stats import fisher_exact, ks_2samp, anderson_ksamp
from numpy.random import multinomial as sample_multinomial
from numpy.random import binomial as sample_binomial


# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
#sample_utils.parse_subject_sample_map

subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sys.stderr.write("Done!\n")

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = ["Bacteroides_vulgatus_57955"]

num_passed_species = 0

for species_name in good_species_list:

    sys.stderr.write("\nProcessing %s...\n" % species_name)

    # First we have to enumerate QP pairs in each cohort
    sys.stderr.write("Enumerating QP pairs...\n")

    # all samples
    all_samples = sample_order_map.keys()

    # list of samples that meet coverage criteria for this species
    highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
    highcoverage_samples = set([highcoverage_sample.decode("utf-8") for highcoverage_sample in highcoverage_samples])

    # list of samples that meet QP criteria for this species
    haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name))
    haploid_samples = set([haploid_sample.decode("utf-8") for haploid_sample in haploid_samples])

    #print len(all_samples), len(highcoverage_samples), len(haploid_samples)

    if len(haploid_samples) < min_haploid_sample_size:
        continue

    all_samples = list(haploid_samples)

    #all_samples = [sample.decode("utf-8") for sample in all_samples]

    same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, all_samples, within_host_type=within_host_type)

    hmp_sample_size = 0

    qp_sample_sets = {cohort: set() for cohort in cohorts}
    qp_counts = {cohort:[0,0,0,0] for cohort in cohorts}

    for sample_pair_idx in range(0,len(same_subject_idxs[0])):

        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]

        sample_i = all_samples[i]
        sample_j = all_samples[j]

        country = sample_country_map[sample_i]

        if country not in countries:
            continue

        # Figure out cohort
        cohort = country_cohort_map[country]

        # Figure out QP status of pair

        if not ((sample_i in highcoverage_samples) and (sample_j in highcoverage_samples)):
            # Both are not highcoverage samples

            if ((sample_i in highcoverage_samples) or (sample_j in highcoverage_samples)):
                # One sample is high coverage
                qp_counts[cohort][0] += 1
            else:
                # Neither sample is high coverage, ignore
                pass

        else:

            # Both are highcoverage samples

            if (sample_i in haploid_samples) and (sample_j in haploid_samples):

                # Both are QP samples!

                qp_counts[cohort][1] += 1
                qp_sample_sets[cohort].add(sample_i)
                qp_sample_sets[cohort].add(sample_j)
                #print sample_i, sample_j

            elif (sample_i not in haploid_samples) and (sample_j not in haploid_samples):
                # pair that is non-QP at both timepoints
                qp_counts[cohort][2] += 1

            else:
                # pair that is QP at one timepoint and non-QP at another
                qp_counts[cohort][3] += 1


    sys.stderr.write("Done!\n")

    #sys.stderr.write("%s\n"%species_name)
    print(species_name)
    for cohort in cohorts:
        print(("%s:" % cohort), qp_counts[cohort][1], "temporal QP pairs.")

    combined_sample_set = set()
    for cohort in cohorts:
        combined_sample_set.update(qp_sample_sets[cohort])
    combined_samples = list(sorted(combined_sample_set))
    combined_sample_idx_map = {combined_samples[i] : i for i in range(0,len(combined_samples))}
    qp_sample_lists = {cohort: list(sorted(qp_sample_sets[cohort])) for cohort in cohorts}

    sample_size = len(qp_sample_sets['hmp'])

    if sample_size < min_sample_size:
        continue

    hmp_species_qp_counts[species_name] = qp_counts['hmp']
    #twin_species_qp_counts[species_name] = qp_counts['twins']

    for cohort in cohorts:
        species_snp_change_distribution[cohort][species_name] = []
        species_snp_nerrs[cohort][species_name] = []
        species_gene_change_distribution[cohort][species_name] = []
        species_gene_nerrs[cohort][species_name] = []

    sys.stderr.write("Proceeding with %d HMP longitudinal comparisons!\n" % (sample_size))

    import calculate_private_snvs
    private_snv_map = calculate_private_snvs.load_private_snv_map(species_name)

    import calculate_snp_prevalences
    snv_freq_map = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)
    snv_freq_keys = list(snv_freq_map.keys())
    snv_freq_values = list(snv_freq_map.values())

    import core_gene_utils
    gene_freq_map = core_gene_utils.parse_gene_freqs(species_name)
    gene_freq_values = numpy.array(list(gene_freq_map.values()))
    gene_freq_weights = gene_freq_values*1.0/gene_freq_values.sum()

    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    import calculate_substitution_rates
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating SNV matrix...\n")
    dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=combined_samples)

    print(snp_mut_opportunity_matrix)


    snp_difference_matrix = snp_mut_difference_matrix+snp_rev_difference_matrix
    snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix
    snp_substitution_rate =     snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    sys.stderr.write("Done!\n")


    sys.stderr.write("Loading gene matrix...\n")
    gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=combined_samples)
    gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
    gene_opportunity_matrix = gene_loss_opportunity_matrix
    gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
    sys.stderr.write("Done!\n")


    print(snv_freq_map[('NC_009614', 5119033)])

    desired_samples = qp_sample_lists['hmp']

    # create presence absence numpy array for private SNVs of each variant type
    for variant_type in ['1D']:

        private_snv_pairwise_sites = []

        for key, snv_dict in private_snv_map.items():

            null_snv_prevalence = [0]* len(desired_samples)

            #print(, snv_dict)

            #if len(snv_dict) == 2:
            #    print(snv_dict)

        #for snv_occupancy in snv_dict:
        #    print(snv_occupancy)



    continue

    ### Need to ask for clarification about opportunity and difference matrices





    sys.stderr.write("Loading pre-computed temporal changes...\n")
    import calculate_temporal_changes
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")

    print(temporal_change_map)

    ### Now loop over different cohorts
    for cohort in cohorts:

        modification_difference_threshold = modification_difference_thresholds[cohort]
        replacement_difference_threshold = replacement_difference_thresholds[cohort]



        same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type)

        #apply_sample_index_map_to_indices(sample_idx_map, idxs):
        #new_idxs = (numpy.array([sample_idx_map[i] for i in idxs[0]]), numpy.array([sample_idx_map[i] for i in idxs[1]]))

        for sample_pair_idx in range(0,len(same_subject_idxs[0])):


            sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]]
            sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]

            i = combined_sample_idx_map[sample_i]
            j = combined_sample_idx_map[sample_j]

            good_idxs = sample_utils.calculate_samples_in_different_subjects( subject_sample_map, combined_samples, sample_i)

            good_idxs *= ( (snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5) )

            if good_idxs.sum() < 1:
                continue

            L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
            nerr = L*perr

            num_mutations = len(mutations)
            num_reversions = len(reversions)
            num_snp_changes = num_mutations+num_reversions

            gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) #, min_normal_copynum = 0.6, max_normal_copynum = 1.2)
            gene_nerr = gene_L*gene_perr
            num_gains = len(gains)
            num_losses = len(losses)
            num_gene_changes = num_gains+num_losses

            if (perr<-0.5) or (gene_perr < -0.5):
                continue

            if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
                continue # Only take things with low-ish FPR

            # Species specific distributions
            species_snp_change_distribution[cohort][species_name].append( num_snp_changes)
            species_snp_nerrs[cohort][species_name].append( nerr)
            species_gene_change_distribution[cohort][species_name].append( num_gene_changes)

            species_gene_nerrs[cohort][species_name].append( gene_nerr)

            # Pooled distributions
            pooled_snp_change_distribution[cohort].append(num_snp_changes)
            pooled_gene_change_distribution[cohort].append(num_gene_changes)

            # Matched between-host samples
            # typical
            pooled_between_snp_change_distribution[cohort].append( choice( snp_difference_matrix[i, good_idxs] ) )
            pooled_between_gene_change_distribution[cohort].append( choice(gene_difference_matrix[i, good_idxs]) )
            # minimum
            pooled_min_between_snp_change_distribution[cohort].append( snp_difference_matrix[i, good_idxs].min() )
            pooled_min_between_gene_change_distribution[cohort].append( gene_difference_matrix[i, good_idxs].min() )

            #if (cohort=='young_twins'):

            #    output_strs[cohort].append("%s, %s: n_snv=%d, n_gene=%d" % (species_name, sample_i, num_snp_changes, num_gene_changes))

            # Store sample names for replacement to see if many species are
            # replaced in the same individual
            if (num_snp_changes>=replacement_difference_threshold):
                sample_pair = (sample_i, sample_j)
                if sample_pair not in replacement_map[cohort]:
                    replacement_map[cohort][sample_pair] = []
                replacement_map[cohort][sample_pair].append(species_name)


            # If deemed a modification, investigate properties of SNVs and genes
            if (num_snp_changes<=modification_difference_threshold):

                if cohort=='twins':
                    pass
                    #output_strs[cohort].append("%s, %s: n_snv=%d, n_gene=%d" % (species_name, sample_i, num_snp_changes, num_gene_changes))

                for snp_change in (mutations+reversions):

                    variant_type = snp_change[3]

                    f = get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map)

                    f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()
                    if variant_type in variant_types:
                        total_freq_snps[cohort][variant_type][f_idx] += 1

                        # Calculate null version based on # of opportunities
                        total_opportunities = 0.0
                        for other_variant_type in variant_types:
                            total_opportunities = opportunity_matrices[other_variant_type][i,j]

                        for other_variant_type in variant_types:
                            total_null_freq_snps[cohort][other_variant_type][f_idx] += opportunity_matrices[other_variant_type][i,j]/total_opportunities

                    total_freq_all_snps[cohort][f_idx] += 1

                    if (sample_i, sample_j, species_name) not in snv_prevalence_count_map[cohort]:
                        snv_prevalence_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_all_snps[cohort])
                        snv_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []

                    snv_prevalence_count_map[cohort][(sample_i, sample_j, species_name)][f_idx] += 1

                    snv_prevalence_map[cohort][(sample_i, sample_j, species_name)].append(f)

                    if variant_type not in variant_type_prevalence_map[cohort]:
                        variant_type_prevalence_map[cohort][variant_type] = []

                    variant_type_prevalence_map[cohort][variant_type].append(f)

                    # Now draw a null prevalence from the genome
                    L = snp_opportunity_matrix[i,j]
                    L_snv = len(snv_freq_map) # A slight overestimate
                    snv_fraction = L_snv*1.0/L
                    num_bootstraps = 10
                    for bootstrap_idx in range(0,num_bootstraps):

                        if random()<snv_fraction:
                            # A polymorphic site

                            random_snv_idx = randint(0,len(snv_freq_keys))
                            random_snv_location = snv_freq_keys[random_snv_idx]
                            f = snv_freq_values[random_snv_idx]

                            rev_f = 1-f

                            if random_snv_location in private_snv_map:
                                # A private SNV. Use private bins
                                # use private
                                f_idx = 0
                                rev_f_idx = -1
                            else:

                                f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()

                                rev_f_idx = ((rev_f>derived_freq_bins[:-1])*(rev_f<=derived_freq_bins[1:])).argmax()


                            # Now add in probability weight
                            total_null_freq_all_snps[cohort][f_idx] += (1-f)*1.0/num_bootstraps
                            total_null_freq_all_snps[cohort][rev_f_idx] += f*1.0/num_bootstraps

                        else:
                            # A truly invariant site
                            total_null_freq_all_snps[cohort][0] += 1.0/num_bootstraps


                for gene_change in gains:
                    gene_name = gene_change[0]

                    if gene_name in gene_freq_map:
                        f = gene_freq_map[gene_name]
                    else:
                        f = 0

                    f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()
                    total_freq_gains[cohort][f_idx] += 1

                    if (sample_i, sample_j, species_name) not in gene_gain_count_map[cohort]:

                        gene_gain_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_gains[cohort])
                        gene_loss_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_losses[cohort])
                        gene_gain_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []
                        gene_loss_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []


                    gene_gain_count_map[cohort][(sample_i, sample_j, species_name)][f_idx] += 1
                    gene_gain_prevalence_map[cohort][(sample_i, sample_j, species_name)].append(f)



                for gene_change in losses:
                    gene_name = gene_change[0]

                    if gene_name in gene_freq_map:
                        f = gene_freq_map[gene_name]
                    else:
                        f = 0

                    f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()
                    total_freq_losses[cohort][f_idx] += 1

                    if (sample_i, sample_j, species_name) not in gene_gain_count_map[cohort]:

                        gene_gain_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_gains[cohort])
                        gene_loss_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_losses[cohort])
                        gene_gain_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []
                        gene_loss_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []

                    gene_loss_count_map[cohort][(sample_i, sample_j, species_name)][f_idx] += 1

                    gene_loss_prevalence_map[cohort][(sample_i, sample_j, species_name)].append(f)


                    num_bootstraps = 10
                    fs = choice(gene_freq_values, size=num_bootstraps, p=gene_freq_weights)
                    for f in fs:
                        f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()
                        total_null_freq_losses[cohort][f_idx] += 1.0/num_bootstraps







#print(variant_type_prevalence_map)

#print(species_snp_change_distribution)

#print( numpy.where(numpy.any((gene_presence_matrix.sum(axis=1) >= 20 ) & (gene_presence_matrix.sum(axis=1) <= len(gene_samples)-20 )  , axis=1)))

# calculate mutual information for all ~1600^2 gene pairs
# not quick, but a lot more optimized than going through all n^2 pairs
# might be as good as it gets

pseudocount = 1 / sum(gene_presence_matrix.sum(axis=0))
normalization_constant = 1 / (len(gene_samples) * (1+pseudocount))
normalization_constant_pairs =  1 / (len(gene_samples) * ((1+pseudocount) ** 2)  )

pseudocount_gene_presence_matrix = gene_presence_matrix + pseudocount
occupancy_probability = normalization_constant * pseudocount_gene_presence_matrix.sum(axis=1)

joint_probability_1_1 = normalization_constant_pairs * numpy.matmul(pseudocount_gene_presence_matrix, pseudocount_gene_presence_matrix.transpose())
joint_probability_1_0 = normalization_constant_pairs * numpy.matmul(pseudocount_gene_presence_matrix, 1+normalization_constant-pseudocount_gene_presence_matrix.transpose())
joint_probability_0_1 = normalization_constant_pairs * numpy.matmul( 1+normalization_constant-pseudocount_gene_presence_matrix, pseudocount_gene_presence_matrix.transpose())
joint_probability_0_0 = normalization_constant_pairs * numpy.matmul( 1+normalization_constant-pseudocount_gene_presence_matrix, 1+normalization_constant-pseudocount_gene_presence_matrix.transpose())


independent_probability_1_1 = numpy.outer(occupancy_probability,occupancy_probability.transpose())
independent_probability_1_0 = numpy.outer(occupancy_probability,(1-occupancy_probability).transpose())
independent_probability_0_1 = numpy.outer((1-occupancy_probability),occupancy_probability.transpose())
independent_probability_0_0 = numpy.outer((1-occupancy_probability),(1-occupancy_probability).transpose())


mutual_information_matrix = joint_probability_1_1 * numpy.log2(numpy.divide(joint_probability_1_1, independent_probability_1_1))
mutual_information_matrix += joint_probability_1_0 * numpy.log2(numpy.divide(joint_probability_1_0, independent_probability_1_0))
mutual_information_matrix += joint_probability_0_1 * numpy.log2(numpy.divide(joint_probability_0_1, independent_probability_0_1))
mutual_information_matrix += joint_probability_0_0 * numpy.log2(numpy.divide(joint_probability_0_0, independent_probability_0_0))

mutual_information_total = (mutual_information_matrix.sum() - numpy.trace(mutual_information_matrix))/2
