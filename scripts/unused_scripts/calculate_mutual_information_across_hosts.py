
from __future__ import division

import config
import parse_midas_data

import os
import parse_HMP_data

import bz2
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
import sample_utils


numpy.random.seed(123456789)


min_coverage = config.min_median_coverage
min_sample_size = 3
min_haploid_sample_size = 10

min_occurances = 10

variant_types = ['1D','4D']

variant_color = {'1D':'r', '4D':'g'}

within_host_type = 'consecutive' # consecutive timepoints
#within_host_type = 'longest' # longest interval available

data_directory = config.data_directory

iterations = 3


slope_null=2

def calculate_sfs_labelled_matrix(species_name, allowed_variant_types=set(['1D','2D','3D','4D']), desired_samples=None, allowed_genes=None, lower_threshold=0.1,upper_threshold=0.9):

    # return site-by-sample matrix for frequency of site within an individual
    # lazy script, not very efficient

    # First write (filtered) genome-wide coverage distribution
    sfs_file = bz2.BZ2File("%ssnps/%s/within_sample_sfs_labelled.txt.bz2" % (data_directory, species_name),"r")
    sfs_file.readline() # header

    #sfs_map = {}
    site_sample_map = {}
    gene_site_map = {}
    samples = []
    for line in sfs_file:
        items = line.split("\t")
        sample = sample_utils.parse_merged_sample_names([items[0].strip()])[0]

        if (sample not in desired_samples) and (type(desired_samples) != type(None)):
            continue


        variant_type = items[1].strip()
        sfs_items = items[2:]


        if variant_type not in allowed_variant_types:
            continue

        #if sample not in sfs_map:
        #    sfs_map[sample] = {}
        #    samples.append(sample)

        if sample not in site_sample_map:
            site_sample_map[sample] = {}


        for sfs_item in sfs_items:
            subitems = sfs_item.split(",")
            D = int(subitems[0])
            A = int(subitems[1])
            n = int(subitems[2])
            reverse_n = float(subitems[3])

            chromosome = subitems[4]
            location = int(subitems[5])
            gene_name = subitems[6]

            variant_type = subitems[7]

            if (gene_name not in allowed_genes) and (type(allowed_genes) != type(None)):
                continue

            if D<0.5:
                continue

            #if (A,D) not in sfs_map[sample]:
            #    sfs_map[sample][(D,A)] = [0,0.0]

            #sfs_map[sample][(D,A)][0] += n
            #sfs_map[sample][(D,A)][1] += reverse_n

            f = A/D

            if (f>lower_threshold) and (f<upper_threshold):

                if f >0.5:
                    f = 1-f

                if gene_name not in gene_site_map:
                    gene_site_map[gene_name] = []

                gene_site_map[gene_name].append(location)

                #if location not in site_sample_map[sample]
                site_sample_map[sample][location] = f


    site_sample_map = {k: v for k, v in site_sample_map.items() if len(v) != 0 }

    all_sites = [v.keys() for k, v in site_sample_map.items()]
    flatten = lambda l: [item for sublist in l for item in sublist]
    all_sites = flatten(all_sites)
    all_sites = list(set(all_sites))
    all_sites.sort()

    site_samples_list = []
    for sample, site_dict in site_sample_map.items():

        site_freqs = [0] *len(all_sites)

        for site_idx, site in enumerate(all_sites):
            if site in site_dict:
                site_freqs[site_idx] = site_dict[site]

        site_samples_list.append(site_freqs)

    site_samples_matrix = numpy.asarray(site_samples_list)

    site_samples_matrix = numpy.transpose(site_samples_matrix)

    return site_samples_matrix, all_sites, site_sample_map.keys(), gene_site_map


def calculate_mutual_information(presence_absence_matrix):

    # calculate mutual information for all ~1600^2 gene pairs
    # not quick, but a lot more optimized than going through all n^2 pairs
    # might be as good as it gets

    number_variables, number_samples = presence_absence_matrix.shape

    pseudocount = 1 / sum(presence_absence_matrix.sum(axis=0))
    normalization_constant = 1 / (number_samples * (1+pseudocount))
    normalization_constant_pairs =  1 / (number_samples * ((1+pseudocount) ** 2)  )

    pseudocount_presence_absence_matrix = presence_absence_matrix + pseudocount

    occupancy_probability = normalization_constant * pseudocount_presence_absence_matrix.sum(axis=1)

    joint_probability_1_1 = normalization_constant_pairs * numpy.matmul(pseudocount_presence_absence_matrix, pseudocount_presence_absence_matrix.transpose())
    joint_probability_1_0 = normalization_constant_pairs * numpy.matmul(pseudocount_presence_absence_matrix, 1+normalization_constant-pseudocount_presence_absence_matrix.transpose())
    joint_probability_0_1 = normalization_constant_pairs * numpy.matmul( 1+normalization_constant-pseudocount_presence_absence_matrix, pseudocount_presence_absence_matrix.transpose())
    joint_probability_0_0 = normalization_constant_pairs * numpy.matmul( 1+normalization_constant-pseudocount_presence_absence_matrix, 1+normalization_constant-pseudocount_presence_absence_matrix.transpose())

    independent_probability_1_1 = numpy.outer(occupancy_probability,occupancy_probability.transpose())
    independent_probability_1_0 = numpy.outer(occupancy_probability,(1-occupancy_probability).transpose())
    independent_probability_0_1 = numpy.outer((1-occupancy_probability),occupancy_probability.transpose())
    independent_probability_0_0 = numpy.outer((1-occupancy_probability),(1-occupancy_probability).transpose())

    mutual_information_matrix = joint_probability_1_1 * numpy.log2(numpy.divide(joint_probability_1_1, independent_probability_1_1))
    mutual_information_matrix += joint_probability_1_0 * numpy.log2(numpy.divide(joint_probability_1_0, independent_probability_1_0))
    mutual_information_matrix += joint_probability_0_1 * numpy.log2(numpy.divide(joint_probability_0_1, independent_probability_0_1))
    mutual_information_matrix += joint_probability_0_0 * numpy.log2(numpy.divide(joint_probability_0_0, independent_probability_0_0))

    mutual_information_total = (mutual_information_matrix.sum() - numpy.trace(mutual_information_matrix))/2

    return mutual_information_total






if __name__=='__main__':


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("--species", help="Name of specific species to run code on", default="all")

    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species=args.species

    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sample_order_map = parse_HMP_data.parse_sample_order_map()
    sample_country_map = parse_HMP_data.parse_sample_country_map()
    sys.stderr.write("Done!\n")

    good_species_list = parse_midas_data.parse_good_species_list()
    if debug:
        #good_species_list = good_species_list[0:4]
        good_species_list = ['Bacteroides_vulgatus_57955']
        #default_num_bootstraps = 100

    elif species !='all':
        good_species_list = [species]


    for species_name in good_species_list:

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
        #snp_samples = [s.decode("utf-8")  for s in snp_samples]

        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough hosts!\n")
            continue
        else:
            sys.stderr.write("Found %d unique hosts!\n" % len(snp_samples))

        sys.stderr.write("Loading core genes...\n")
        core_genes = parse_midas_data.load_core_genes(species_name)
        sys.stderr.write("Done! %d core genes\n" % len(core_genes))
        allowed_genes = core_genes


        site_samples_matrix_nonsynonymous, sites_nonsynonymous, snp_samples_nonsynonymous, gene_site_map_nonsynonymous = calculate_sfs_labelled_matrix(species_name, allowed_variant_types=set(['1D']), desired_samples=snp_samples, allowed_genes=allowed_genes)
        site_samples_matrix_synonymous, sites_synonymous, snp_samples_synonymous, gene_site_map_synonymous = calculate_sfs_labelled_matrix(species_name, allowed_variant_types=set(['4D']), desired_samples=snp_samples, allowed_genes=allowed_genes)

        # presence absence matrix
        #site_samples_matrix_presence_abscence_nonsynonymous = numpy.where(site_samples_matrix_nonsynonymous > 0, 1, 0)
        #site_samples_matrix_presence_abscence_synonymous = numpy.where(site_samples_matrix_synonymous > 0, 1, 0)

        #site_counts = site_samples_matrix_presence_abscence.sum(axis=1)
        # remove sites that occur less than min_occurances times
        #site_samples_matrix_presence_abscence = site_samples_matrix_presence_abscence[site_counts[site_counts>=min_occurances],:]
        # remove samples

        #prevalence_nonsynonymous = site_samples_matrix_presence_abscence_nonsynonymous.sum(axis=1) / site_samples_matrix_presence_abscence_nonsynonymous.shape[1]
        #prevalence_synonymous = site_samples_matrix_presence_abscence_synonymous.sum(axis=1) / site_samples_matrix_presence_abscence_synonymous.shape[1]

        # quick and dirty sfs FFD, taylors law

        import matplotlib.pyplot as plt


        fig, ax = plt.subplots(figsize=(4,4))

        #ax.plot([0.05,0.5],[0.05,0.5], lw=3,ls='--',c='k',zorder=1)

        for variant_type in variant_types:

            site_samples_matrix, sites, snp_samples, gene_site_map = calculate_sfs_labelled_matrix(species_name, allowed_variant_types=set([variant_type]), desired_samples=snp_samples, allowed_genes=allowed_genes)

            #site_samples_matrix_presence_abscence = numpy.where(site_samples_matrix > 0, 1, 0)

            #prevalence = site_samples_matrix_presence_abscence.sum(axis=1) / site_samples_matrix_presence_abscence.shape[1]

            means = []
            variances = []

            prevalences = []

            for i in range(site_samples_matrix.shape[0]):
                frequencies = site_samples_matrix[i,:][site_samples_matrix[i,:] > 0]

                if len(frequencies) > 3:
                    #print(sites[i] , frequencies)
                    means.append(numpy.mean(frequencies))
                    variances.append(numpy.var(frequencies))

                    prevalences.append(len(frequencies) /  site_samples_matrix.shape[1])


            ax.scatter(means, prevalences, alpha=0.05, c=variant_color[variant_type], label=variant_type)#, c='#87CEEB')


        ax.set_xlabel('Average SNP MAF', fontsize=12)
        ax.set_ylabel('SNP prevalence', fontsize=10)

        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)

        ax.legend(loc="lower right", fontsize=8)

        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        fig.savefig(config.analysis_directory + "prevalence_mean_frequency.png", format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
        plt.close()





        continue


        cutoffs = [0,0.12,1]

        for variant_type in variant_types:

            site_samples_matrix, sites, snp_samples, gene_site_map = calculate_sfs_labelled_matrix(species_name, allowed_variant_types=set([variant_type]), desired_samples=snp_samples, allowed_genes=allowed_genes)
            site_samples_matrix_presence_abscence = numpy.where(site_samples_matrix > 0, 1, 0)

            prevalence = site_samples_matrix_presence_abscence.sum(axis=1) / site_samples_matrix_presence_abscence.shape[1]

            for cutoff_idx in range(len(cutoffs)-1):
                #prevalence_cutoff_idx =


                site_samples_matrix_presence_abscence_cutoff = site_samples_matrix_presence_abscence[(prevalence >= cutoffs[cutoff_idx]) & (prevalence < cutoffs[cutoff_idx+1]),:]

                mutual_information_total = calculate_mutual_information(site_samples_matrix_presence_abscence_cutoff)

                mutual_information_total_null_list = []

                for i in range(iterations):

                    print(i)

                    site_samples_matrix_presence_abscence_copy = site_samples_matrix_presence_abscence_cutoff.copy()

                    for j in range(site_samples_matrix_presence_abscence_copy.shape[0]):
                        numpy.random.shuffle(site_samples_matrix_presence_abscence_copy[j,:])

                    mutual_information_total_null = calculate_mutual_information(site_samples_matrix_presence_abscence_copy)

                    mutual_information_total_null_list.append(mutual_information_total_null)

                mutual_information_total_null_list.sort()

                standardized_mutual_information = (mutual_information_total - numpy.mean(mutual_information_total_null_list)) /  numpy.std(mutual_information_total_null_list)

                print(variant_type, cutoffs[cutoff_idx], mutual_information_total, standardized_mutual_information)
