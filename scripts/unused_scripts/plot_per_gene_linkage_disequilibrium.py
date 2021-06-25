import os
import gzip

import scipy
import matplotlib.pyplot as plt
import numpy
import itertools

import config

import parse_midas_data
import calculate_linkage_disequilibria


intermediate_filename_template = '%s%s.txt.gz'

clade_types = ['all','largest_clade']
variant_types = ['4D','1D']


good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = ['Bacteroides_vulgatus_57955']

ld_per_gene_directory = '%slinkage_disequilibria_per_gene/' % (parse_midas_data.data_directory)

per_gene_ld_dict_all_species = {}

for clade_type in clade_types:
    # clade_type ==> species ==> sigma N / sigma S
    per_gene_ld_dict_all_species[clade_type] = {}

for species_name in good_species_list:



    intermediate_filename = intermediate_filename_template % (ld_per_gene_directory, species_name)

    if not os.path.isfile(intermediate_filename):
        continue

    ld_dict = calculate_linkage_disequilibria.load_per_gene_ld_map(species_name)

    per_gene_ld_dict = {}

    print(species_name)


    #for clade_type in clade_types:
    #    for variant_type in variant_types:
    #        per_gene_ld_dict[(clade_type,variant_type)] = {}

    for clade_type in clade_types:

        per_gene_ld_dict_all_species[clade_type][species_name] = []

        for variant_type in variant_types:

            for key, items in ld_dict[(clade_type,variant_type)].items():

                # correct for distance
                # different genes have different sizes ==> can't assume equal sampling of distances across gene
                # multiply LD by (obs - min) / (max - min)
                distances = numpy.asarray(items['distances'])
                rsquared_numerators = numpy.asarray(items['rsquared_numerators'])
                rsquared_denominators = numpy.asarray(items['rsquared_denominators'])

                rsquards = rsquared_numerators/rsquared_denominators

                #remove negative LD values, ask Nandita about this
                # why are there LD values equal to one?
                rsquards_no_negative = rsquards[(rsquards>0) & (distances >= items['minimum_distance']) & (rsquards!=1)]
                distances_no_negative = distances[(rsquards>0) & (distances >= items['minimum_distance']) & (rsquards!=1) ]


                size_corrected_distances = (distances_no_negative - items['minimum_distance'] + 1) / (items['gene_size']-1 - items['minimum_distance'] +1)

                rsquards_no_negative_size_corrected = rsquards_no_negative*size_corrected_distances

                #if len(rsquards_no_negative_size_corrected) < 10:
                #    continue

                if key not in per_gene_ld_dict:
                    per_gene_ld_dict[key] = {}

                per_gene_ld_dict[key][(clade_type,variant_type)] = rsquards_no_negative_size_corrected


    for gene_name, gene_name_dict in per_gene_ld_dict.items():

        for clade_type_i in clade_types:

            if ((clade_type_i, '1D') in gene_name_dict.keys()) and ((clade_type_i, '4D')  in gene_name_dict.keys()):

                rsquards_no_negative_size_corrected_1D = gene_name_dict[(clade_type_i, '1D')]
                rsquards_no_negative_size_corrected_4D = gene_name_dict[(clade_type_i, '4D')]

                if (len(rsquards_no_negative_size_corrected_1D) >= 10) and (len(rsquards_no_negative_size_corrected_4D) >= 10):

                    #per_gene_ld_dict_all_species[clade_type_i][species_name].append(numpy.mean(rsquards_no_negative_size_corrected_1D)/numpy.mean(rsquards_no_negative_size_corrected_4D))

                    per_gene_ld_dict_all_species[clade_type_i][species_name].append(numpy.mean(rsquards_no_negative_size_corrected_1D) - numpy.mean(rsquards_no_negative_size_corrected_4D))


# plot
#fig = plt.figure(figsize = (4, 4)) #

fig, ax = plt.subplots(figsize = (4, 4))

ax.axvline(x=0, color='k', linestyle=':', alpha = 0.8, zorder=2)


for species, rsquard_ratios in  per_gene_ld_dict_all_species['largest_clade'].items():

    rsquard_ratios = numpy.asarray(rsquard_ratios)

    if len(rsquard_ratios) < 100:
        continue

    #low, high = numpy.floor(rsquard_ratios.min()), numpy.ceil(rsquard_ratios.max())
    #bins = numpy.linspace(low, high, high - low + 1)

    #weights = np.ones_like(myarray)/float(len(myarray))

    #plt.bar(edges[:-1], y_values*binWidth, binWidth)

    #results, edges = numpy.histogram(myarray, normed=True)


    y_values, bin_edges = numpy.histogram(rsquard_ratios, bins=20, density=True)

    binWidth = bin_edges[1] - bin_edges[0]

    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])


    #plt.plot(bin_centers, y_values, '-',lw=2, alpha=0.7)#, label=data_set)
    #ax.plot(bin_centers, y_values*binWidth, '-',lw=2, alpha=0.7)#, label=data_set)

    ax.scatter(bin_centers, y_values*binWidth, alpha=0.4)




#wspace=0.3, hspace=0.3
#plt.xscale('log')
#

ax.set_xlabel('$\overline{\sigma}^{2}_{d, N} - \overline{\sigma}^{2}_{d, S}$' , fontsize = 12)

ax.set_ylabel('Fraction of genes' , fontsize = 12)

fig.savefig("%sper_gene_linkage_disequilibrium.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

plt.close()
