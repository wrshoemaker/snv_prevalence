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
import os
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

species_to_run = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                    'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                    'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507',
                    'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318',
                    'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208',
                    'Dialister_invisus_61905', 'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985',
                    'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']

# to-redo....
# Bacteroides_ovatus_58035



def get_mapgd_data(species_name, samples, whitelisted_genes):

    gene_location_name_dict = calculate_predicted_prevalence_mapgd.get_gene_location_name_dict(species_name)

    mapgd_dict = {}

    # parse MAPGD output
    mapgd_samples = []
    chromosome_location_tuples = []
    for sample in samples:

        #mapgd_file_path = '/u/home/w/wrshoema/project-ngarud/HMP_MAPGD/midas_output_bam/%s/%s_sorted.pol' % (species_name, sample)
        mapgd_file_path = '%smidas_output_bam/%s/%s_sorted.pol' % (config.data_directory, species_name, sample)

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
                mapgd_dict[gene_name][location]['location_tuples'] = []

            mapgd_dict[gene_name][location]['samples'].append(sample)
            mapgd_dict[gene_name][location]['major_alleles'].append(major)
            mapgd_dict[gene_name][location]['minor_alleles'].append(minor)
            mapgd_dict[gene_name][location]['frequencies'].append(frequency)
            mapgd_dict[gene_name][location]['coverages'].append(coverage)
            mapgd_dict[gene_name][location]['location_tuples'].append((chromosome, location))

            #if (chromosome, location) not in chromosome_location_tuples:
            chromosome_location_tuples.append((chromosome, location))

        mapgd_file.close()
    chromosome_location_tuples = list(set(chromosome_location_tuples))

    return mapgd_samples, mapgd_dict, chromosome_location_tuples




def make_mapgd_output_dict():

    for species_name in species_to_run:

        sys.stderr.write("%s\n" % species_name)

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

        sys.stderr.write("Loading MAPGD data...\n")
        mapgd_samples, mapgd_dict, chromosome_location_tuples = get_mapgd_data(species_name, samples, core_genes)
        mapgd_samples = numpy.asarray(mapgd_samples)

        mapgd_output_dict = {}
        mapgd_output_dict['frequency_minor_allele'] = []
        mapgd_output_dict['coverage_minor_allele'] = []
        mapgd_output_dict['coverage_total'] = []

        for gene in mapgd_dict.keys():
            for location in mapgd_dict[gene].keys():
                frequencies = mapgd_dict[gene][location]['frequencies']
                coverages = mapgd_dict[gene][location]['coverages']
                coverages = [int(c) for c in coverages]

                for frequency_idx, frequency in enumerate(frequencies):
                    minor_frequency = min([frequency, 1-frequency])
                    minor_coverage = coverages[frequency_idx]*minor_frequency

                    mapgd_output_dict['frequency_minor_allele'].append(minor_frequency)
                    mapgd_output_dict['coverage_minor_allele'].append(minor_coverage)
                    mapgd_output_dict['coverage_total'].append(coverages[frequency_idx])


        # parse through, get alternate coverage counts, minimum frequencies

        intermediate_filename = "%smapgd_output_dicts/%s.dat" % (config.data_directory, species_name)

        sys.stderr.write("Saving mapgd summary dict...\n")
        with open(intermediate_filename, 'wb') as handle:
            pickle.dump(mapgd_output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



# Alistipes_finegoldii_56071

def get_bam_coverage_data():

    bam_species = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                    'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                    'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035',
                    'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']

    bam_species = ['Alistipes_finegoldii_56071']

    for species_name in bam_species:

        sys.stderr.write("%s\n" % species_name)

        sys.stderr.write("Loading whitelisted genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)

        # Holds panel wide prevalence for each species
        #os.system('mkdir -p %ssnp_prevalences' % config.data_directory)
        snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        #if os.path.isfile(snp_file_path) == False:
        #    continue

        # Open post-processed MIDAS output
        snp_file =  bz2.BZ2File(snp_file_path, "r")
        line = snp_file.readline() # header
        items = line.split()[1:]
        snp_file.close()

        samples = numpy.array([item.strip() for item in items])

        sys.stderr.write("Loading MAPGD data...\n")
        mapgd_samples, mapgd_dict, chromosome_location_tuples = get_mapgd_data(species_name, samples, core_genes)

        bam_coverage_directory = '%sbam_coverage/%s/' % (config.data_directory, species_name)

        bam_coverage_dict = {}

        n_sites = 0

        for filename in os.listdir(bam_coverage_directory):

            if filename.endswith("_sorted.gz"):

                filepath = os.path.join(bam_coverage_directory, filename)

                # ignore empty files
                if os.stat(filepath).st_size == 20:
                    continue

                f = gzip.open(filepath, 'r')
                for line in f:

                    if n_sites > 10000:
                        break

                    line = line.strip().split('\t')
                    chromosome_location_tuple = (line[0], int(line[1]))
                    if chromosome_location_tuple in chromosome_location_tuples:
                        coverage = int(line[2])

                        if chromosome_location_tuple not in bam_coverage_dict:
                            bam_coverage_dict[chromosome_location_tuple] = []

                        print(coverage)

                        bam_coverage_dict[chromosome_location_tuple].append(coverage)

                        n_sites += 1

                f.close()


        intermediate_filename = "%sbam_coverage_dicts/%s.dat" % (config.data_directory, species_name)
        sys.stderr.write("Saving bam coverage dict...\n")
        with open(intermediate_filename, 'wb') as handle:
            pickle.dump(bam_coverage_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


get_bam_coverage_data()
