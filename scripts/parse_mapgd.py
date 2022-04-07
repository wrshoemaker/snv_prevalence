import numpy
import sys
import bz2
import gzip
import os
import os.path
import stats_utils
from math import floor, ceil
import gene_diversity_utils

import config
import sample_utils
import parse_midas_data
import parse_HMP_data
import diversity_utils
import calculate_substitution_rates
import clade_utils

import pickle

low_divergence_threshold = config.between_low_divergence_threshold
#low_divergence_threshold = 5e-04 # this was picked by looking at inflection point of dN/dS vs dS plot
min_sample_size = config.between_host_min_sample_size
min_ld_sample_size = config.between_host_ld_min_sample_size
allowed_variant_types = set(['1D','4D'])

#use annotated_snps.txt.bz2 to get annotations

# base code off of this


def make_snp_annotation_map(species_name, initial_line_number=0):


    # Open post-processed MIDAS output
    snp_file =  bz2.BZ2File("%ssnps/%s/annotated_snps.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")

    line = snp_file.readline() # header
    items = line.split()[1:]

    num_sites_processed = 0
    line_number = -1
    final_line_number = -1
    previous_gene_name = ""
    gene_name = ""

    annotation_map = {}

    for line in snp_file:

        line_number += 1
        previous_gene_name = gene_name

        if line_number < initial_line_number:
            continue

        items = line.split()
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = long(info_items[1])
        gene_name = info_items[2]
        variant_type = info_items[3]

        if chromosome not in annotation_map:
            annotation_map[chromosome] = {}

        annotation_map[chromosome][location] = items[0]

    snp_file.close()


    #print(annotation_map)

    # dump annotation map

    output_name = "%ssnps/%s/annotated_snp_map.dat" % (parse_midas_data.data_directory, species_name)

    with open(output_name, 'wb') as handle:
        pickle.dump(annotation_map, handle,  protocol=pickle.HIGHEST_PROTOCOL)


def load_snp_annotation_map(species_name):

    output_name = "%ssnps/%s/annotated_snp_map.dat" % (parse_midas_data.data_directory, species_name)

    with open(output_name, 'rb') as handle:
        annotation_map = pickle.load(handle)

    return annotation_map


def calculate_prevalences_from_mapgd(species_name):

    annotation_map = load_snp_annotation_map(species_name)

    directory_name = "%smapgd_output/%s/" % (parse_midas_data.data_directory, species_name)

    site_map = {}

    samples = []

    for file_name in os.listdir(directory_name):
        if file_name.endswith(".pol") == False:
            continue

        sample = file_name.split('.')[0].split('_')[0]

        samples.append(sample)

        with open("%s%s" % (directory_name, file_name)) as file:

            line_number = 0

            for line in file.readlines():

                line_number += 1

                if (line_number == 1) or (line_number == 2):
                    continue

                #line = file.readline()
                line = line.strip().split()

                if len(line) < 8:
                    continue

                chromosome = line[0]
                location = int(line[1])
                major_allele = line[3]
                minor_allele = line[4]
                coverage = int(line[5])
                frequency = float(line[7].split('/')[0])

                if chromosome not in site_map:
                    site_map[chromosome] = {}

                if location not in site_map[chromosome]:
                    site_map[chromosome][location] = {}
                    #site_map[chromosome][location]['samples'] = []
                    #site_map[chromosome][location]['frequencies'] = []

                    #site_map[chromosome][location]['major_alleles'] = []
                    #site_map[chromosome][location]['minor_alleles'] = []


                #site_map[chromosome][location]['samples'].append(sample)
                #site_map[chromosome][location]['frequencies'].append(frequency)
                #site_map[chromosome][location]['major_alleles'].append(major_allele)
                #site_map[chromosome][location]['minor_alleles'].append(minor_allele)

                site_map[chromosome][location][sample] = {}
                site_map[chromosome][location][sample]['frequency'] = frequency
                site_map[chromosome][location][sample]['major_allele'] = major_allele
                site_map[chromosome][location][sample]['minor_allele'] = minor_allele
                site_map[chromosome][location][sample]['coverage'] = coverage


    output_file = bz2.BZ2File("%ssnps/%s/annotated_snps_mapgd.txt.bz2" % (parse_midas_data.data_directory, species_name),"w")

    header_line = samples[:]
    header_line.insert(0,'site_id')

    output_file.write("\t".join(header_line))



    for chromosome, chromosome_dict in site_map.iteritems():

        if chromosome not in annotation_map:
            continue


        for site, site_dict in chromosome_dict.iteritems():

            if site not in annotation_map[chromosome]:
                continue

            annotation  = annotation_map[chromosome][site]


            frequencies_sites = [ str(site_dict[sample]['frequency']) if sample in site_dict.keys() else str(float(0)) for sample in samples ]

            frequencies_sites.insert(0,annotation)

            #frequencies_sites_join =
            output_file.write("\t".join(frequencies_sites))
            output_file.write("\n")


    output_file.close()
    sys.stderr.write("MAPGD annotated file done!\n")




















###############################################################################
#
# Loads list of SNPs and counts of target sites from annotated SNPs file
#
# returns (lots of things, see below)
#
###############################################################################
def parse_snps(species_name, debug=False, allowed_samples=[], allowed_genes=[], allowed_variant_types=['1D','2D','3D','4D'], initial_line_number=0, chunk_size=1000000000):

    import calculate_snp_prevalences
    # Load population freqs (for polarization purposes)
    population_freqs = calculate_snp_prevalences.parse_population_freqs(species_name, polarize_by_consensus=False)

    # Open post-processed MIDAS output
    #snp_file =  bz2.BZ2File("%ssnps/%s/annotated_snps.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")
    snp_file =  bz2.BZ2File("%ssnps/%s/annotated_snps_mapgd.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")

    line = snp_file.readline() # header
    items = line.split()[1:]
    samples = sample_utils.parse_merged_sample_names(items)

    if len(allowed_samples)==0:
        allowed_sample_set = set(samples)
    else:
        allowed_sample_set = (set(allowed_samples) & set(samples))

    allowed_genes = set(allowed_genes)
    allowed_variant_types = set(allowed_variant_types)

    # This is a hack because there were some mistaken repeats in an old data file
    # should be able to remove later
    seen_samples = set()
    desired_sample_idxs = []
    for sample in allowed_samples:
        if (sample in allowed_sample_set) and (sample not in seen_samples):
            desired_sample_idxs.append( numpy.nonzero(samples==sample)[0][0] )
        else:
            pass

        seen_samples.add(sample)

    desired_sample_idxs = numpy.array(desired_sample_idxs)

    desired_samples = samples[desired_sample_idxs]

    # map from gene_name -> var_type -> (list of locations, matrix of allele counts)
    allele_counts_map = {}
    # map from gene_name -> var_type -> (location, sample x sample matrix of whether both samples can be called at that site)
    passed_sites_map = {}

    num_sites_processed = 0
    line_number = -1
    final_line_number = -1
    previous_gene_name = ""
    gene_name = ""
    for line in snp_file:

        line_number += 1
        previous_gene_name = gene_name

        if line_number < initial_line_number:
            continue

        items = line.split()
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = long(info_items[1])
        gene_name = info_items[2]
        variant_type = info_items[3]

        if len(info_items) > 5: # for backwards compatability
            polarization = info_items[4]
            pvalue = float(info_items[5])
        else:
            polarization="R" # not correct, but avoids a crash
            pvalue = float(info_items[4])


        if num_sites_processed >= chunk_size and gene_name!=previous_gene_name:
            # We are done for now!
            final_line_number = line_number
            break


        if not variant_type in allowed_variant_types:
            continue

        if len(allowed_genes)>0 and (not gene_name in allowed_genes):
            continue

        # Load alt and depth counts
        # Load alt and depth counts
        #alts = []
        #depths = []

        frequencies = []


        for idx in desired_sample_idxs:
            item = items[1+idx]
            frequencies.append(float(item))
            #subitems = item.split(",")
            #alts.append(float(subitems[0]))
            #depths.append(float(subitems[1]))
        #alts = numpy.array(alts)
        #depths = numpy.array(depths)

        frequencies = numpy.asarray(frequencies)


        # polarize
        if (chromosome, location) in population_freqs:
            population_freq = population_freqs[(chromosome, location)]
        else:
            population_freq = 0

        # polarize SFS according to population freq
        if population_freq>0.5:
            #alts = depths-alts
            frequencies = 1-frequencies
            polarization = 'A'

        #passed_sites = (depths>0)*1.0
        passed_sites = numpy.asarray([True]*len(frequencies))*1.0
        if gene_name not in passed_sites_map:
            passed_sites_map[gene_name] = {v: {'location': (chromosome,location), 'sites': numpy.zeros((len(desired_samples), len(desired_samples)))} for v in allowed_variant_types}

            #allele_counts_map[gene_name] = {v: {'locations':[], 'alleles':[]} for v in allowed_variant_types}
            allele_counts_map[gene_name] = {v: {'locations':[], 'frequencies':[]} for v in allowed_variant_types}

        passed_sites_map[gene_name][variant_type]['sites'] += passed_sites[:,None]*passed_sites[None,:]

        # zero out non-passed sites
        # (shouldn't be needed anymore)
        #alts = alts*passed_sites
        #depths = depths*passed_sites

        #frequencies = frequencies*frequencies


        # calculate whether SNP has passed
        #alt_threshold = numpy.ceil(depths*config.parse_snps_min_freq)+0.5 #at least one read above 5%.
        #snp_passed = ((alts>alt_threshold).sum()>0) and (pvalue<0.05)
        snp_passed = ((frequencies>0).sum()>0)


        if snp_passed:
            #allele_counts = numpy.transpose(numpy.array([alts,depths-alts]))

            allele_counts_map[gene_name][variant_type]['locations'].append((chromosome, location))
            #allele_counts_map[gene_name][variant_type]['alleles'].append(allele_counts)
            allele_counts_map[gene_name][variant_type]['frequencies'].append(frequencies)

            num_sites_processed+=1

            if num_sites_processed>0 and num_sites_processed%1000==0:
                sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))
                if debug:
                    break

    snp_file.close()


    for gene_name in passed_sites_map.keys():
        for variant_type in passed_sites_map[gene_name].keys():

            #allele_counts_map[gene_name][variant_type]['alleles'] = numpy.array(allele_counts_map[gene_name][variant_type]['alleles'])
            allele_counts_map[gene_name][variant_type]['frequencies'] = numpy.array(allele_counts_map[gene_name][variant_type]['frequencies'])

    return desired_samples, allele_counts_map, passed_sites_map, final_line_number




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




    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")


    #sys.exit("stop script")


    good_species_list = parse_midas_data.parse_good_species_list()
    if species!='all':
        good_species_list = [species]
    if debug and len(good_species_list)>3.5:
        #good_species_list = good_species_list[:3]
        good_species_list = ['Bacteroides_vulgatus_57955']


    for species_name in good_species_list:

        #make_snp_annotation_map(species_name)

        #calculate_prevalences_from_mapgd(species_name)

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


        final_line_number = 0
        while final_line_number >= 0:

            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_snps(species_name, debug=debug,         allowed_variant_types=allowed_variant_types, allowed_samples=snp_samples,allowed_genes=core_genes, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))

            print(allele_counts_map)
