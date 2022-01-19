
import sys
import numpy
import config

import diversity_utils
import parse_midas_data
import parse_HMP_data
import calculate_substitution_rates
import sample_utils

#import matplotlib.pyplot as plt


from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster, to_tree

intermediate_filename_template = '%s%s%s.txt'

low_divergence_threshold = config.between_low_divergence_threshold
#low_divergence_threshold = 5e-04 # this was picked by looking at inflection point of dN/dS vs dS plot
min_sample_size = config.between_host_min_sample_size
min_ld_sample_size = config.between_host_ld_min_sample_size

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick



#good_species_list = ['Escherichia_coli_58110', 'Eubacterium_rectale_56927']

good_species_list = parse_midas_data.parse_good_species_list()
#good_species_list = ['Eubacterium_rectale_56927']
#species_name = 'Eubacterium_rectale_56927'
for species_name in good_species_list:
    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sample_continent_map = parse_HMP_data.parse_sample_continent_map()

    sys.stderr.write("Loading haploid samples...\n")
        # Only plot samples above a certain depth threshold that are "haploids"
    snp_samples = diversity_utils.calculate_haploid_samples(species_name)

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
    snp_samples = dummy_samples

    snp_substitution_matrix = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))

    snp_substitution_rate = snp_substitution_matrix
    snp_substitution_rate = numpy.clip(snp_substitution_rate,1e-11,10)

    sys.stderr.write("Calculating UPGMA dendrogram...\n")
    # calculate compressed distance matrix suitable for agglomerative clustering
    Y = []
    for i in xrange(0,snp_substitution_rate.shape[0]):
        for j in xrange(i+1,snp_substitution_rate.shape[1]):
            Y.append(snp_substitution_rate[i,j])
    Y = numpy.array(Y)
    Z = linkage(Y, method='average')
    c, coph_dists = cophenet(Z, Y)
    #ddata = dendrogram(Z, no_plot=True)
    sys.stderr.write("Done! cophenetic correlation: %g\n" % c)

    tree = to_tree(Z,False)
    sys.stderr.write("Writing Newick file...\n")

    newick = getNewick(tree, "", tree.dist, dummy_samples)

    intermediate_filename = intermediate_filename_template % (config.data_directory, 'cladogram_newick/', species_name)
    output_file = open(intermediate_filename,"w")

    output_file.write(newick)

    output_file.close()

    sys.stderr.write("Done!\n")
