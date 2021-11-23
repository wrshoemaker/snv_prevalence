import sys, os
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria

import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#import calculate_snv_distances
import figure_utils
from math import log10,ceil
#import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.random import randint, multinomial
import parse_HMP_data

import figure_utils
import ld_utils

import matplotlib as mpl

import calculate_coprevalence


bin_width_exponent = 0.3
variant_types = ['4D','1D']

color_dict = {'4D': '#87CEEB', '1D':'#FF6347'}


sys.stderr.write("Loading sample metadata...\n")
#subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species_list = parse_midas_data.parse_good_species_list()



for species_name in good_species_list:

    coprevalence_directory = '%scoprevalence/%s.txt.gz' % (parse_midas_data.data_directory, species_name)

    if os.path.isfile(coprevalence_directory) == False:
        continue

    coprevalence_dict = calculate_coprevalence.load_ld_map(species_name)

    fig, ax = plt.subplots(figsize=(4,4))

    for variant_type in variant_types:

        coprevalence_dict_4D = coprevalence_dict[('all', variant_type)]

        distances = coprevalence_dict_4D[0]

        intragene_rsquared_numerators = coprevalence_dict_4D[1]

        intragene_rsquared_denominators = coprevalence_dict_4D[2]

        intragene_rsquared_ratio = intragene_rsquared_numerators/intragene_rsquared_denominators


        distances_to_plot = distances[(intragene_rsquared_ratio>0) & numpy.logical_not(numpy.isnan(intragene_rsquared_ratio))]
        intragene_rsquared_ratio_to_plot = intragene_rsquared_ratio[(intragene_rsquared_ratio>0) & numpy.logical_not(numpy.isnan(intragene_rsquared_ratio))]

        #distances_to_plot = distances[intragene_rsquared_denominators>0]
        #intragene_rsquared_numerators_to_plot = intragene_rsquared_numerators[intragene_rsquared_denominators > 0]
        #intragene_rsquared_denominators_to_plot = intragene_rsquared_denominators[intragene_rsquared_denominators > 0]

        #x = x[]


        color = color_dict[variant_type]

        #print(intragene_rsquared_ratio_to_plot)

        #ax.plot(distances_to_plot, intragene_rsquared_numerators_to_plot/intragene_rsquared_denominators_to_plot, color=color, alpha = 0.7, label = variant_type)
        ax.plot(distances_to_plot, intragene_rsquared_ratio_to_plot, color=color, alpha = 0.7, label = variant_type)


    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

    ax.set_xlabel('Distance between SNPs', fontsize=12)
    ax.set_ylabel('Coprevalence', fontsize=12)

    ax.legend(loc="lower left")


    ax.set_title(species_name)

    fig.tight_layout()
    fig.savefig("%scoprevalence_vs_distance/%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



    #load_ld_map(good_species_list[0])
