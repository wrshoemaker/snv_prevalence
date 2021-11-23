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
from numpy.random import randint, choice, multinomial


import midas_db_utils

import parse_mapgd

import parse_HMP_data

intermediate_filename_template = '%s%s_f0_%f_k_%d.txt.gz'



low_divergence_threshold = config.between_low_divergence_threshold
#low_divergence_threshold = 5e-04 # this was picked by looking at inflection point of dN/dS vs dS plot
min_sample_size = config.between_host_min_sample_size
min_ld_sample_size = config.between_host_ld_min_sample_size
allowed_variant_types = set(['1D','4D'])


#f0_range = numpy.logspace(numpy.log10(0.01), numpy.log10(0.3), num=30, endpoint=True, base=10.0)
#f0_range = numpy.logspace(numpy.log10(0.01), numpy.log10(0.3), num=2, endpoint=True, base=10.0)



def load_ld_map(species_name, f0, ld_moment):

    ld_map = {}

    ld_directory = '%scoprevalence_f0/' % (parse_midas_data.data_directory)

    #intermediate_filename = intermediate_filename_template % (ld_directory, species_name)
    intermediate_filename = intermediate_filename_template % (ld_directory, species_name, f0, ld_moment)


    if not os.path.isfile(intermediate_filename):
        return ld_map

    file = gzip.open(intermediate_filename,"r")
    header_line = file.readline() # header
    header_items = header_line.split(",")

    distance_strs = [item.split(":")[-1] for item in header_items[4:]]

    distances = []
    intragene_idxs = []

    intergene_distances = []
    intergene_idxs = []

    control_idx = -1

    for i in xrange(0,len(distance_strs)-1):

        if distance_strs[i].startswith('g'):
            # an intergene distance

            intergene_idxs.append(i)
            intergene_distances.append(long(distance_strs[i][1:]))
        else:
            # an intragene distance
            intragene_idxs.append(i)
            distances.append(float(distance_strs[i]))

    distances = numpy.array(distances)
    intragene_idxs = numpy.array(intragene_idxs)

    intergene_distances = numpy.array(intergene_distances)
    intergene_idxs = numpy.array(intergene_idxs)

    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue

        clade_type = items[1].strip()
        variant_type = items[2].strip()
        pi = float(items[3])

        rsquared_numerators = []
        rsquared_denominators = []
        #lds = []
        counts = []
        for item in items[4:]:
            subitems = item.split(":")
            rsquared_numerators.append(float(subitems[0]))
            rsquared_denominators.append(float(subitems[1]))
            counts.append(float(subitems[2]))

        rsquared_numerators = numpy.array(rsquared_numerators)
        rsquared_denominators = numpy.array(rsquared_denominators)
        counts = numpy.array(counts)

        #rsquared_numerators = rsquared_numerators[rsquared_denominators!=0]
        #rsquared_denominators = rsquared_denominators[rsquared_denominators!=0]
        #counts = counts[rsquared_denominators!=0]

        #lds = rsquared_numerators/rsquared_denominators

        control_numerator = rsquared_numerators[control_idx]
        control_denominator = rsquared_denominators[control_idx]
        control_count = counts[control_idx]

        control_ld = control_numerator/control_denominator

        intragene_rsquared_numerators = rsquared_numerators[intragene_idxs]
        intragene_rsquared_denominators = rsquared_denominators[intragene_idxs]
        intragene_counts = counts[intragene_idxs]

        intergene_rsquared_numerators = rsquared_numerators[intergene_idxs]
        intergene_rsquared_denominators = rsquared_denominators[intergene_idxs]
        intergene_counts = counts[intergene_idxs]

        ld_map[(clade_type, variant_type)] = (distances, intragene_rsquared_numerators, intragene_rsquared_denominators, intragene_counts, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_counts, control_numerator, control_denominator, control_count, pi)

    file.close()

    return ld_map







def calculate_ld_dict(good_species_list, subject_sample_map, f0, ld_moment, bin_width_exponent=0.1, variant_types = ['4D','1D'], debug=False):

    ld_dict_all_species = {}
    passed_species = []
    sample_sizes = []

    for species_name in good_species_list:

        ld_directory = '%scoprevalence_f0/' % (parse_midas_data.data_directory)

        #intermediate_filename_template = '%s%s.txt.gz'
        #intermediate_filename = intermediate_filename_template % (ld_directory, species_name)
        intermediate_filename_template = '%s%s_f0_%f_k_%d.txt.gz'
        intermediate_filename = intermediate_filename_template % (ld_directory, species_name, f0, ld_moment)


        if os.path.exists(intermediate_filename) == False:
            continue

        sys.stderr.write("Loading haploid samples...\n")
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

        sys.stderr.write("Calculating unique samples...\n")

        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


        # Load precomputed LD
        ld_map = load_ld_map(species_name, f0, ld_moment)
        ld_dict = {}

        for variant_type in variant_types:

            ld_dict[variant_type] = {}

            if len(ld_map)>0:

                # fix this when you rerun the code
                #distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerators, control_rsquared_denominators, control_n, pi = ld_map[('largest_clade',variant_type)]
                distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerators, control_rsquared_denominators, control_n, pi = ld_map[('all',variant_type)]

                if True:
                    passed_species.append(species_name)
                    sample_sizes.append(len(snp_samples))
                else:
                    sys.stderr.write("%s intergene LD too high: %g (%g)\n" % (species_name, control_rsquared, rsquareds[0]))


            all_distances, all_rsquared_numerators, all_rsquared_denominators, all_ns, all_intergene_distances, all_intergene_rsquared_numerators, all_intergene_rsquared_denominators, all_intergene_ns, all_control_rsquared_numerator, all_control_rsquared_denominator, all_control_n, all_pi = ld_map[('all',variant_type)]
            all_control_rsquared = all_control_rsquared_numerator/all_control_rsquared_denominator

            distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerator, control_rsquared_denominator, control_n, pi = ld_map[('largest_clade',variant_type)]
            control_rsquared = control_rsquared_numerator/control_rsquared_denominator

            #distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerator, control_rsquared_denominator, control_n, pi = ld_map[('all',variant_type)]
            #control_rsquared = control_rsquared_numerator/control_rsquared_denominator

            # smooth this stuff:
            smoothed_distances = distances
            #window_width = 10**(0.1)
            window_width = 10**(bin_width_exponent)

            dmins = smoothed_distances/(window_width**0.5)
            dmaxs = smoothed_distances*(window_width**0.5)

            smoothed_rsquared_numerators = []
            smoothed_rsquared_denominators = []
            smoothed_counts = []

            all_smoothed_rsquared_numerators = []
            all_smoothed_rsquared_denominators = []
            all_smoothed_counts = []

            for dmin,dmax in zip(dmins,dmaxs):
                binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
                binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
                binned_counts = ns[(distances>=dmin)*(distances<=dmax)]

                smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
                smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
                smoothed_counts.append( binned_counts.sum() )

                binned_numerators = all_rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
                binned_denominators = all_rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
                binned_counts = all_ns[(distances>=dmin)*(distances<=dmax)]
                all_smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
                all_smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
                all_smoothed_counts.append( binned_counts.sum() )


            smoothed_rsquared_numerators = numpy.array( smoothed_rsquared_numerators )
            smoothed_rsquared_denominators = numpy.array( smoothed_rsquared_denominators )
            smoothed_counts = numpy.array( smoothed_counts )

            all_smoothed_rsquared_numerators = numpy.array( all_smoothed_rsquared_numerators )
            all_smoothed_rsquared_denominators = numpy.array( all_smoothed_rsquared_denominators )
            all_smoothed_counts = numpy.array( all_smoothed_counts )

            early_distances = distances[distances<101]
            early_rsquareds = rsquared_numerators[distances<101]*1.0/rsquared_denominators[distances<101]
            early_ns = ns[distances<101]

            early_distances = early_distances[early_ns>0.5]
            early_rsquareds = early_rsquareds[early_ns>0.5]
            early_ns = early_ns[early_ns>0.5]


            distances = smoothed_distances
            rsquareds = smoothed_rsquared_numerators/(smoothed_rsquared_denominators)
            ns = smoothed_counts
            distances = distances[ns>0]
            rsquareds = rsquareds[ns>0]


            # adding this code below to fix error after filtering "distances"
            rsquared_numerators = rsquared_numerators[ns>0]
            rsquared_denominators = rsquared_denominators[ns>0]

            # then filter ns

            ns = ns[ns>0]

            all_distances = smoothed_distances
            #all_distances = dmins
            all_rsquareds = all_smoothed_rsquared_numerators/(all_smoothed_rsquared_denominators)
            all_ns = all_smoothed_counts

            all_distances = all_distances[all_ns>0]
            all_rsquareds = all_rsquareds[all_ns>0]
            all_ns = all_ns[all_ns>0]

            num_bootstraps = 10

            bootstrapped_sigmasquareds = [] # will eventually be a matrix where first index is window_idx and second index is bootstrap index (sorted from lowest to highest)
            # Estimate bootstrap intervals for focal species only
            for dmin,dmax in zip(dmins,dmaxs):
                binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
                binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
                binned_counts = ns[(distances>=dmin)*(distances<=dmax)]

                total_pairs = binned_counts.sum()

                upper_rsquareds = []
                lower_rsquareds = []

                if total_pairs>0:

                    if True: #len(binned_counts)>1:
                        #print total_pairs
                        #print binned_counts
                        ps = binned_counts*1.0/total_pairs

                        window_bootstrapped_countss = multinomial(total_pairs,ps,size=num_bootstraps)

                        #print window_bootstrapped_countss.shape
                        window_bootstrapped_numerators = (window_bootstrapped_countss*binned_numerators[None,:]).sum(axis=1)*1.0/total_pairs
                        window_bootstrapped_denominators = (window_bootstrapped_countss*binned_denominators[None,:]).sum(axis=1)*1.0/total_pairs

                        window_bootstrapped_sigmasquareds = window_bootstrapped_numerators/window_bootstrapped_denominators

                        #print window_bootstrapped_sigmasquareds.shape
                        window_bootstrapped_sigmasquareds.sort()

                        bootstrapped_sigmasquareds.append(window_bootstrapped_sigmasquareds)

                        #print total_pairs

                    else:
                        bootstrapped_sigmasquareds.append([binned_numerators/binned_denominators]*num_bootstraps)

                else:

                    bootstrapped_sigmasquareds.append([-1]*num_bootstraps)


            upper_rsquareds = numpy.array([bootstrapped_sigmasquareds[window_idx][int(num_bootstraps*0.95)] for window_idx in range(0,len(bootstrapped_sigmasquareds))])
            lower_rsquareds = numpy.array([bootstrapped_sigmasquareds[window_idx][int(num_bootstraps*0.05)] for window_idx in range(0,len(bootstrapped_sigmasquareds))])

            #print upper_rsquareds-lower_rsquareds

            good_distances = (upper_rsquareds>=-0.5)*(lower_rsquareds>=-0.5)

            #theory_ls = numpy.logspace(0,log10(distances[-1]),100)
            #theory_NRs = theory_ls/200.0
            #theory_rsquareds = (10+2*theory_NRs)/(22+26*theory_NRs+4*theory_NRs*theory_NRs)
            ld_dict[variant_type]['all_distances'] = all_distances
            ld_dict[variant_type]['distances'] = distances
            ld_dict[variant_type]['good_distances'] = good_distances
            ld_dict[variant_type]['early_distances'] = early_distances

            ld_dict[variant_type]['all_control_rsquared'] = all_control_rsquared

            ld_dict[variant_type]['rsquareds'] = rsquareds
            ld_dict[variant_type]['all_rsquareds'] = all_rsquareds
            ld_dict[variant_type]['lower_rsquareds'] = lower_rsquareds
            ld_dict[variant_type]['upper_rsquareds'] = upper_rsquareds
            ld_dict[variant_type]['early_rsquareds'] = early_rsquareds
            ld_dict[variant_type]['control_rsquared'] = control_rsquared


        ld_dict_all_species[species_name] = ld_dict


    return ld_dict_all_species




def get_f0_range(n_range=40, largest_clade=False):

    # Load subject and sample metadata
    #sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()

    #sys.stderr.write("Done!\n")

    good_species_list_all = parse_midas_data.parse_good_species_list()
    #if species!='all':
    #    good_species_list = [species]
    #if debug and len(good_species_list)>3.5:
    #    #good_species_list = good_species_list[:3]
    #    good_species_list = ['Bacteroides_vulgatus_57955']

    n_samples_list = []

    for species_name in good_species_list_all:

        #sys.stderr.write("Loading haploid samples...\n")
        # Only plot samples above a certain depth threshold that are "haploids"
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

        if len(snp_samples) < min_sample_size:
            #sys.stderr.write("Not enough haploid samples!\n")
            continue
        #else:
        #    sys.stderr.write("Found %d haploid samples!\n" % len(snp_samples))

        #sys.stderr.write("Calculating unique hosts...\n")
        # Only consider one sample per person
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]

        if len(snp_samples) < min_sample_size:
            #sys.stderr.write("Not enough hosts!\n")
            continue
        #else:
        #    sys.stderr.write("Found %d unique hosts!\n" % len(snp_samples))

        # Load divergence matrices
        #sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        #sys.stderr.write("Calculating matrix...\n")
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = numpy.array(dummy_samples)
        substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        #sys.stderr.write("Done!\n")

        #sys.stderr.write("Clustering samples with low divergence...\n")
        coarse_grained_idxs, coarse_grained_cluster_list = clade_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold) # NRG: what is this returning?

        coarse_grained_samples = snp_samples[coarse_grained_idxs]
        clade_sets = clade_utils.load_manual_clades(species_name)

        #sys.stderr.write("%d samples remaining after clustering!\n" % len(coarse_grained_samples))

        if len(clade_sets)==0:
            continue

        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)

        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])

        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]

        largest_clade_set = set(largest_clade_samples)

        #sys.stderr.write("Top level clades: %d clades, sizes: %s\n" % (len(clade_sets), str(clade_sizes)))
        #sys.stderr.write("Max clade size: %d\n" % len(largest_clade_samples))

        snp_samples = coarse_grained_samples

        if len(largest_clade_samples) < min_ld_sample_size:
            #sys.stderr.write("Not enough ld samples!\n")
            continue
        #else:
        #    sys.stderr.write("Proceeding with %d coarse-grained samples in largest clade!\n" % len(largest_clade_samples))

        if largest_clade == False:
            n_samples_list.append(len(snp_samples))

        else:
            n_samples_list.append(len(largest_clade_samples))


    f0_range = numpy.logspace(numpy.log10(1/min(n_samples_list)), numpy.log10(0.1), num=n_range, endpoint=True, base=10.0)

    return f0_range



if __name__=='__main__':


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("--species", help="Name of specific species to run code on", default="all")
    parser.add_argument("--ld_moment", type=int, help="Moment of LD to estimate", default=2)

    #parser.add_argument('--condition_on_freq', help='Conditions on allele frequency', action='store_true')
    #parser.add_argument('--low', help='allele frequency lower bound', default=0.0, type=float)
    #parser.add_argument('--high', help='allele frequency upper bound', default=1.0, type=float)

    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species=args.species
    ld_moment=args.ld_moment



    ld_directory = '%scoprevalence_f0/' % (parse_midas_data.data_directory)


    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")

    good_species_list = parse_midas_data.parse_good_species_list()
    if species!='all':
        good_species_list = [species]
    if debug and len(good_species_list)>3.5:
        #good_species_list = good_species_list[:3]
        good_species_list = ['Bacteroides_vulgatus_57955']



    #good_species_list=['Bacteroides_vulgatus_57955']
    # better binning scheme (multiple of 3)
    distance_bin_locations = numpy.arange(0,1002)*3.0
    distance_bins = numpy.arange(-1,1002)*3+1.5
    distance_bins[0] = 0 # no such thing as negative distance
    distance_bins[1] = 2.5 # want at least one codon separation
    distance_bins[-1] = 1e09 # catch everything

    neighbor_distances = numpy.array([1,2,3,4,5])


    distance_strs = ["LD_N:LD_D:%g" % d for d in distance_bin_locations[1:-1]] # N=numerator and D=denominator
    distance_strs = distance_strs+["LD_N:LD_D:g%d" % nd for nd in neighbor_distances]+["LD_N:LD_D:intergene"]

    # header of the output file.
    record_strs = [", ".join(['Species', 'CladeType', 'VariantType', 'Pi']+distance_strs)]


    os.system('mkdir -p %s' % ld_directory)

    #f0_range = get_f0_range()

    f0_range = [0.05]


    for species_name in good_species_list:

        sys.stderr.write("Loading haploid samples...\n")
        # Only plot samples above a certain depth threshold that are "haploids"
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        else:
            sys.stderr.write("Found %d haploid samples!\n" % len(snp_samples))

        sys.stderr.write("Calculating unique hosts...\n")
        # Only consider one sample per person
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough hosts!\n")
            continue
        else:
            sys.stderr.write("Found %d unique hosts!\n" % len(snp_samples))

        # Load divergence matrices
        sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        sys.stderr.write("Calculating matrix...\n")
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = numpy.array(dummy_samples)
        substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        sys.stderr.write("Done!\n")

        sys.stderr.write("Clustering samples with low divergence...\n")
        coarse_grained_idxs, coarse_grained_cluster_list = clade_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold) # NRG: what is this returning?

        coarse_grained_samples = snp_samples[coarse_grained_idxs]
        clade_sets = clade_utils.load_manual_clades(species_name)

        sys.stderr.write("%d samples remaining after clustering!\n" % len(coarse_grained_samples))

        if len(clade_sets)==0:
            continue

        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)

        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])

        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]

        largest_clade_set = set(largest_clade_samples)

        sys.stderr.write("Top level clades: %d clades, sizes: %s\n" % (len(clade_sets), str(clade_sizes)))
        sys.stderr.write("Max clade size: %d\n" % len(largest_clade_samples))

        snp_samples = coarse_grained_samples

        if len(largest_clade_samples) < min_ld_sample_size:
            sys.stderr.write("Not enough ld samples!\n")
            continue
        else:
            sys.stderr.write("Proceeding with %d coarse-grained samples in largest clade!\n" % len(largest_clade_samples))




        # Analyze SNPs, looping over chunk sizes.
        # Clunky, but necessary to limit memory usage on cluster

        sys.stderr.write("Loading core genes...\n")
        core_genes = parse_midas_data.load_core_genes(species_name)
        sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        sys.stderr.write("(core genes only...)\n")

        #f0_range =  numpy.logspace(numpy.log10(1/len(largest_clade_samples)), numpy.log10( (len(largest_clade_samples)-1) / len(largest_clade_samples) ), num=20, base=10, endpoint=True)

        clade_types = ['all','largest_clade']
        #clade_types = ['all']
        variant_types = ['4D','1D']

        for f0 in f0_range:

            # add to these nested dictionaries

            binned_rsquared_numerators = {}
            binned_rsquared_denominators = {}
            binned_counts = {}

            neighboring_gene_rsquared_numerators = {}
            neighboring_gene_rsquared_denominators = {}
            neighboring_gene_counts = {}

            # total_control=between genes.
            total_control_rsquared_numerators = {}
            total_control_rsquared_denominators = {}
            total_control_counts = {}

            total_control_genes = {}



            gene_size_dict = midas_db_utils.get_reference_gene_sizes(species_name)


            for clade_type in clade_types:
                for variant_type in variant_types:

                    binned_rsquared_numerators[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)
                    binned_rsquared_denominators[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)
                    binned_counts[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)

                    neighboring_gene_rsquared_numerators[(clade_type,variant_type)] = numpy.zeros_like(neighbor_distances)*1.0
                    neighboring_gene_rsquared_denominators[ (clade_type,variant_type)] = numpy.zeros_like(neighbor_distances)*1.0
                    neighboring_gene_counts[(clade_type,variant_type)] = numpy.zeros_like(neighbor_distances)*1.0

                    total_control_rsquared_numerators[(clade_type,variant_type)] = 0
                    total_control_rsquared_denominators[(clade_type,variant_type)] = 0
                    total_control_counts[(clade_type,variant_type)] = 0

                    total_control_genes[(clade_type,variant_type)] = []



            final_line_number = 0
            while final_line_number >= 0:

                sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
                #snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug,         allowed_variant_types=allowed_variant_types, allowed_samples=snp_samples,allowed_genes=core_genes, chunk_size=chunk_size,initial_line_number=final_line_number)
                snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_mapgd.parse_snps(species_name, debug=debug,         allowed_variant_types=allowed_variant_types, allowed_samples=snp_samples,allowed_genes=core_genes, chunk_size=chunk_size,initial_line_number=final_line_number)
                sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))

                largest_clade_idxs = numpy.array([sample in largest_clade_set for sample in snp_samples])

                sys.stderr.write("Calculating LD...\n")
                for clade_type in clade_types:

                    for variant_type in variant_types:

                        for gene_name in allele_counts_map.keys():

                            if gene_name not in core_genes:
                                continue

                            locations = numpy.array([location for chromosome, location in allele_counts_map[gene_name][variant_type]['locations']])*1.0
                            #allele_counts = allele_counts_map[gene_name][variant_type]['alleles']
                            allele_counts = allele_counts_map[gene_name][variant_type]['frequencies']


                            if len(allele_counts)==0:
                                # no diversity to look at!
                                continue

                            target_chromosome = allele_counts_map[gene_name][variant_type]['locations'][0][0]

                            if clade_type=='largest_clade':
                                # Now restrict to largest clade
                                #allele_counts = allele_counts[:,largest_clade_idxs,:]
                                allele_counts = allele_counts[largest_clade_idxs]


                            #compute the distances between all pairs of sites
                            # None in the two index positions results in a transpose of the vector relative to each other
                            # Subtraction between the two vectors results in pairwise subtraction of each element in each vector.
                            distances = numpy.fabs(locations[:,None]-locations[None,:])

                            print(allele_counts)



                            #rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k(allele_counts, allele_counts, f0=f0, k=ld_moment)
                            #rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k_prevalence(allele_counts, allele_counts, f0, k=ld_moment)
                            rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k_prevalence_mapgd(allele_counts, allele_counts, f0, k=ld_moment)

                            print("what")
                            print(rsquared_numerators, rsquared_denominators)

                            neighbor_rsquared_numeratorss = [[] for d in neighbor_distances]
                            neighbor_rsquared_denominatorss = [[] for d in neighbor_distances]

                            for neighbor_distance_idx in xrange(0,len(neighbor_distances)):

                                neighbor_distance = neighbor_distances[neighbor_distance_idx]

                                gene_name_items = gene_name.split(".")
                                gene_peg_number = long(gene_name.split(".")[-1])
                                nearest_gene_peg_numbers = [gene_peg_number-neighbor_distance,gene_peg_number+neighbor_distance]
                                neighboring_genes = [".".join(gene_name_items[:-1]+[str(n)]) for n in nearest_gene_peg_numbers]

                                for neighboring_gene_name in neighboring_genes:

                                    # first make sure it's a real gene
                                    if neighboring_gene_name not in allele_counts_map:
                                        continue

                                    if neighboring_gene_name not in core_genes:
                                        continue

                                    #neighboring_allele_counts = allele_counts_map[neighboring_gene_name][variant_type]['alleles']
                                    neighboring_allele_counts = allele_counts_map[neighboring_gene_name][variant_type]['frequencies']

                                    # then make sure it has some variants
                                    if len(neighboring_allele_counts)==0:
                                        continue

                                    neighboring_chromosome = allele_counts_map[neighboring_gene_name][variant_type]['locations'][0][0]

                                    if neighboring_chromosome!=target_chromosome:
                                        continue

                                    if clade_type=='largest_clade':
                                    # Now restrict to largest clade
                                        #neighboring_allele_counts = neighboring_allele_counts[:,largest_clade_idxs,:]
                                        neighboring_allele_counts = neighboring_allele_counts[largest_clade_idxs]


                                    #chunk_rsquared_numerators, chunk_rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k(allele_counts, neighboring_allele_counts, f0=f0, k=ld_moment)
                                    chunk_rsquared_numerators, chunk_rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k_prevalence_mapgd(allele_counts, neighboring_allele_counts, f0, k=ld_moment)

                                    neighbor_rsquared_numeratorss[ neighbor_distance_idx].extend( chunk_rsquared_numerators.flatten() )
                                    neighbor_rsquared_denominatorss[ neighbor_distance_idx].extend( chunk_rsquared_denominators.flatten() )

                                neighbor_rsquared_numeratorss[ neighbor_distance_idx] = numpy.array( neighbor_rsquared_numeratorss[neighbor_distance_idx] )

                                neighbor_rsquared_denominatorss[ neighbor_distance_idx] = numpy.array( neighbor_rsquared_denominatorss[neighbor_distance_idx] )

                            # pick a random gene somewhere else as a control
                            # 10 to 1 control to regular
                            control_rsquared_numerators = []
                            control_rsquared_denominators = []
                            gene_peg_number = long(gene_name.split(".")[-1])

                            for control_idx in xrange(0,10):

                                control_gene_name = gene_name
                                control_allele_counts = []

                                # get the right gene name
                                while True:

                                    control_gene_name = choice(allele_counts_map.keys())

                                    if control_gene_name not in core_genes:
                                        continue


                                    control_gene_peg_number = long(control_gene_name.split(".")[-1])

                                    control_allele_counts = allele_counts_map[control_gene_name][variant_type]['frequencies']

                                    if len(control_allele_counts)==0:
                                        continue

                                    # make sure you don't have one too close by!
                                    if (fabs(control_gene_peg_number - gene_peg_number) < 5.5):
                                        continue

                                    if clade_type=='largest_clade':
                                    # Now restrict to largest clade
                                        #control_allele_counts = control_allele_counts[:,largest_clade_idxs,:]
                                        control_allele_counts = control_allele_counts[largest_clade_idxs]

                                    break


                                #control_gene_rsquared_numerators, control_gene_rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k(allele_counts, control_allele_counts, f0=f0, k=ld_moment)
                                #control_gene_rsquared_numerators, control_gene_rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k_prevalence(allele_counts, control_allele_counts, f0, k=ld_moment)
                                control_gene_rsquared_numerators, control_gene_rsquared_denominators = diversity_utils.calculate_unbiased_sigma_k_prevalence_mapgd(allele_counts, control_allele_counts, f0, k=ld_moment)

                                control_rsquared_numerators.extend( control_gene_rsquared_numerators.flatten() )
                                control_rsquared_denominators.extend( control_gene_rsquared_denominators.flatten() )

                            control_rsquared_numerators = numpy.array( control_rsquared_numerators )
                            control_rsquared_denominators = numpy.array( control_rsquared_denominators )

                            # get the indices of the upper diagonal of the distance matrix
                            # numpy triu_indices returns upper diagnonal including diagonal
                            # the 1 inside the function excludes diagonal. Diagnonal has distance of zero.
                            print(distances)
                            desired_idxs = numpy.triu_indices(distances.shape[0],1)

                            print(desired_idxs)


                            # fetch the distances and rsquared vals corresponding to the upper diagonal.
                            distances = distances[desired_idxs]
                            print(desired_idxs)
                            print(rsquared_numerators)
                            rsquared_numerators = rsquared_numerators[desired_idxs]
                            rsquared_denominators = rsquared_denominators[desired_idxs]

                            print(rsquared_numerators)

                            # fetch entries where denominator != 0 (remember, denominator=pa*(1-pa)*pb*(1-pb). If zero, then at least one site is invariant)


                            #distances = distances[rsquared_denominators>1e-09]
                            #rsquared_numerators = rsquared_numerators[rsquared_denominators>1e-09]
                            #rsquared_denominators = rsquared_denominators[rsquared_denominators>1e-09]

                            distances = distances[((rsquared_denominators>1e-09) & ( ~numpy.isnan(rsquared_denominators)) & (numpy.isfinite(rsquared_denominators)))]
                            rsquared_numerators = rsquared_numerators[((rsquared_denominators>1e-09) & ( ~numpy.isnan(rsquared_denominators)) & (numpy.isfinite(rsquared_denominators)))]
                            rsquared_denominators = rsquared_denominators[((rsquared_denominators>1e-09) & ( ~numpy.isnan(rsquared_denominators)) & (numpy.isfinite(rsquared_denominators)))]




                            if len(distances) == 0:
                                continue



                            # numpy.digitize: For each distance value, return the bin index it belongs to in distances_bins.
                            bin_idxs = numpy.digitize(distances,bins=distance_bins)-1


                            for i in xrange(0,len(bin_idxs)):

                                binned_counts[(clade_type,variant_type)][bin_idxs[i]] += 1
                                binned_rsquared_numerators[(clade_type,variant_type)][bin_idxs[i]] += rsquared_numerators[i]
                                binned_rsquared_denominators[(clade_type,variant_type)][bin_idxs[i]] += rsquared_denominators[i]


                                #print(binned_rsquared_denominators[(clade_type,variant_type)][bin_idxs[i]])

                            for i in xrange(0,len(neighbor_distances)):
                                good_idxs = ( (neighbor_rsquared_denominatorss[i]>1e-09) & ( ~numpy.isnan(neighbor_rsquared_denominatorss[i])) & (numpy.isfinite(neighbor_rsquared_denominatorss[i])) )

                                neighboring_gene_counts[(clade_type,variant_type)][i] += good_idxs.sum()

                                neighboring_gene_rsquared_numerators[ (clade_type,variant_type)][i] += neighbor_rsquared_numeratorss[i][good_idxs].sum()

                                neighboring_gene_rsquared_denominators[ (clade_type,variant_type)][i] += neighbor_rsquared_denominatorss[i][good_idxs].sum()


                            #total_control_counts[(clade_type,variant_type)] += (control_rsquared_denominators>1e-09).sum()
                            #total_control_rsquared_numerators[(clade_type,variant_type)] += control_rsquared_numerators[control_rsquared_denominators>1e-09].sum()
                            #total_control_rsquared_denominators[(clade_type,variant_type)] += control_rsquared_denominators[control_rsquared_denominators>1e-09].sum()

                            total_control_counts[(clade_type,variant_type)] += ((control_rsquared_denominators>1e-09) & ( ~numpy.isnan(control_rsquared_denominators)) & (numpy.isfinite(control_rsquared_denominators))).sum()
                            total_control_rsquared_numerators[(clade_type,variant_type)] += control_rsquared_numerators[((control_rsquared_denominators>1e-09) & ( ~numpy.isnan(control_rsquared_denominators)) & (numpy.isfinite(control_rsquared_denominators)))].sum()
                            total_control_rsquared_denominators[(clade_type,variant_type)] += control_rsquared_denominators[((control_rsquared_denominators>1e-09) & ( ~numpy.isnan(control_rsquared_denominators)) & (numpy.isfinite(control_rsquared_denominators)))].sum()




                            total_control_genes[(clade_type,variant_type)].append(gene_name)



            for clade_type in clade_types:
                for variant_type in variant_types:

                    desired_samples = snp_samples
                    if clade_type=='largest_clade':
                        desired_samples = largest_clade_samples

                    # Calculate pi!
                    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, variant_type, allowed_samples=desired_samples)

                    substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))

                    iu = numpy.triu_indices(substitution_rate.shape[0], 1)

                    pi = numpy.median(substitution_rate[iu])


                    binned_rsquareds = binned_rsquared_numerators[(clade_type,variant_type)]*1.0/(binned_rsquared_denominators[(clade_type,variant_type)] + (binned_rsquared_denominators[(clade_type,variant_type)] == 0))

                    control_rsquareds = total_control_rsquared_numerators[(clade_type,variant_type)]*1.0/(total_control_rsquared_denominators[(clade_type,variant_type)]+(total_control_rsquared_denominators[(clade_type,variant_type)]==0))

                    rsquared_strs = ["%g:%g:%d" % (rsquared_numerator, rsquared_denominator, count) for rsquared_numerator, rsquared_denominator, count in zip(binned_rsquared_numerators[(clade_type,variant_type)], binned_rsquared_denominators[(clade_type,variant_type)], binned_counts[(clade_type,variant_type)])[1:-1]]

                    gene_rsquared_strs = ["%g:%g:%d" % (rsquared_numerator, rsquared_denominator, count) for rsquared_numerator, rsquared_denominator, count in zip(neighboring_gene_rsquared_numerators[(clade_type,variant_type)], neighboring_gene_rsquared_denominators[(clade_type,variant_type)], neighboring_gene_counts[(clade_type,variant_type)])]

                    control_rsquared_str = "%g:%g:%d" % (total_control_rsquared_numerators[(clade_type,variant_type)], total_control_rsquared_denominators[(clade_type,variant_type)], total_control_counts[(clade_type,variant_type)])

                    pi_str = str(pi)


                    record_str = ", ".join([species_name, clade_type, variant_type, pi_str]+rsquared_strs+gene_rsquared_strs+[control_rsquared_str])

                    record_strs.append(record_str)

            sys.stderr.write("Done with %s!\n" % species_name)


            sys.stderr.write("Writing intermediate file...\n")
            intermediate_filename = intermediate_filename_template % (ld_directory, species_name, f0, ld_moment)
            file = gzip.open(intermediate_filename,"w")
            record_str = "\n".join(record_strs)
            file.write(record_str)
            file.close()
            sys.stderr.write("Done!\n")




            sys.stderr.write("Done looping over species!\n")

            sys.stderr.write("Testing loading...\n")
            ld_map = load_ld_map(good_species_list[0], f0_range[0], ld_moment)

            sys.stderr.write("Done!\n")
