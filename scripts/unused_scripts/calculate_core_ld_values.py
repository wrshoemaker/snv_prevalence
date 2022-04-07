import sys, os
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria

#import calculate_snv_distances
from math import log10,ceil
from numpy.random import randint, multinomial
import parse_HMP_data

focal_speciess = ['Bacteroides_vulgatus_57955', 'Akkermansia_muciniphila_55290']

species_name='Bacteroides_vulgatus_57955'

focal_colors = ['b','g']

variant_types = ['4D','1D']

color_dict = {'4D':'b', '1D':'r'}

passed_species = []
sample_sizes = []


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize



debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species_list = parse_midas_data.parse_good_species_list()




ld_dict_all_species = {}

output_file = open('/u/home/w/wrshoema/project-ngarud/genome_wide_ld.txt', "w")

output_file.write("\t".join(['Species', 'Variant type', 'Clade cutoff type', 'Genome wide LD']))
output_file.write("\n")
for species_name in good_species_list:

    ld_directory = '%slinkage_disequilibria/' % (parse_midas_data.data_directory)
    intermediate_filename_template = '%s%s.txt.gz'
    intermediate_filename = intermediate_filename_template % (ld_directory, species_name)

    if os.path.exists(intermediate_filename) == False:
        continue

    sys.stderr.write("Loading haploid samples...\n")
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

    sys.stderr.write("Calculating unique samples...\n")

    snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


    # Load precomputed LD
    ld_map = calculate_linkage_disequilibria.load_ld_map(species_name)
    ld_dict = {}

    for variant_type in variant_types:

        ld_dict[variant_type] = {}


        if len(ld_map)>0:

            distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerators, control_rsquared_denominators, control_n, pi = ld_map[('largest_clade',variant_type)]

            if True:
                passed_species.append(species_name)
                sample_sizes.append(len(snp_samples))
            else:
                sys.stderr.write("%s intergene LD too high: %g (%g)\n" % (species_name, control_rsquared, rsquareds[0]))


        all_distances, all_rsquared_numerators, all_rsquared_denominators, all_ns, all_intergene_distances, all_intergene_rsquared_numerators, all_intergene_rsquared_denominators, all_intergene_ns, all_control_rsquared_numerator, all_control_rsquared_denominator, all_control_n, all_pi = ld_map[('all',variant_type)]
        all_control_rsquared = all_control_rsquared_numerator/all_control_rsquared_denominator

        distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerator, control_rsquared_denominator, control_n, pi = ld_map[('largest_clade',variant_type)]
        control_rsquared = control_rsquared_numerator/control_rsquared_denominator


        #ld_dict[variant_type]['all_control_rsquared'] = all_control_rsquared
        #ld_dict[variant_type]['control_rsquared'] = control_rsquared

        output_file.write("\t".join([species_name, variant_type, 'all', str(all_control_rsquared)]))
        output_file.write("\n")
        output_file.write("\t".join([species_name, variant_type, 'largest_clade', str(control_rsquared)]))
        output_file.write("\n")



output_file.close()
