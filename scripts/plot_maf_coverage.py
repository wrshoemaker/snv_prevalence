from __future__ import division
import sys
import os
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
import plot_utils
import prevalence_utils

from scipy.stats import gamma
import scipy.special
from scipy.integrate import quad


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

import calculate_predicted_prevalence_mapgd


species_color_map, ordered_species_list = plot_utils.get_species_color_map()


species_to_run = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                    'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                    'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035',
                    'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318',
                    'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208',
                    'Dialister_invisus_61905', 'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985',
                    'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']



max_cov_minor = 20
max_cov_total = 40
range_minor = list(range(max_cov_minor+1))
range_total = list(range(max_cov_total+1))


#fig, ax = plt.subplots(figsize=(4.5,4))
gs = gridspec.GridSpec(nrows=1, ncols=3)
fig = plt.figure(figsize = (12, 4))

ax_total_coverage = fig.add_subplot(gs[0, 0])
ax_minor_coverage = fig.add_subplot(gs[0, 1])
ax_scatter = fig.add_subplot(gs[0, 2])


ax_total_coverage.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_total_coverage.transAxes)
ax_minor_coverage.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_minor_coverage.transAxes)
ax_scatter.text(-0.1, 1.04, prevalence_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_scatter.transAxes)


coverage_total_counts_all = []
coverage_minor_counts_all = []

#species_to_run = [species_to_run[0]]

for species_name in species_to_run:

    print(species_name)

    intermediate_filename_template = config.data_directory+"mapgd_output_dicts/%s.dat"
    intermediate_filename = intermediate_filename_template % species_name

    if os.path.isfile(intermediate_filename) == False:
        continue

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)

    frequency_minor_allele = b['frequency_minor_allele']
    coverage_minor_allele = b['coverage_minor_allele']
    coverage_total = b['coverage_total']

    coverage_minor_allele = [int(round(c)) for c in coverage_minor_allele]

    frequency_minor_allele = numpy.asarray(frequency_minor_allele)
    coverage_minor_allele = numpy.asarray(coverage_minor_allele)
    coverage_total = numpy.asarray(coverage_total)


    #frequency_minor_allele = numpy.log10(frequency_minor_allele)

    #ax_freq.hist(frequency_minor_allele, bins=30, density=True, histtype='step')
    #ax_coverage.hist(coverage_minor_allele, bins=30, density=True, histtype='step')

    coverage_total_counts = numpy.asarray([sum(coverage_total==i) for i in range(max_cov_total+1)])
    coverage_total_counts = coverage_total_counts/len(coverage_total)

    coverage_minor_allele_counts = numpy.asarray([sum(coverage_minor_allele==i) for i in range(max_cov_minor+1)])
    coverage_minor_allele_counts = coverage_minor_allele_counts/len(coverage_minor_allele)


    ax_total_coverage.scatter(range_total, coverage_total_counts, color=species_color_map[species_name], s=20, alpha=0.8)
    ax_minor_coverage.scatter(range_minor, coverage_minor_allele_counts, color=species_color_map[species_name], s=20, alpha=0.8)


    coverage_total_to_plot = coverage_total[coverage_total<=max_cov_total]


    range_scatter = range(20, 100+1)
    #for i in range_scatter:

    #    print(i, coverage_minor_allele[coverage_total==i])
    coverage_minor_allele_mean = numpy.asarray([numpy.median(coverage_minor_allele[coverage_total==i]) for i in range_scatter])
    #print(coverage_minor_allele_mean)

    ax_scatter.scatter(range_scatter, coverage_minor_allele_mean, color=species_color_map[species_name], s=20, alpha=0.8)


    coverage_total_counts_all.extend(coverage_total_counts.tolist())
    coverage_minor_counts_all.extend(coverage_minor_allele_counts.tolist())


ax_total_coverage.set_xlim(1 , max_cov_total)
ax_minor_coverage.set_xlim(1 , max_cov_minor)
#ax_scatter.set_xlim(1 , max_cov_minor)

#ax.set_ylim(0, 0.18)

ax_total_coverage.set_ylim(0, max(coverage_total_counts_all)*1.1)
ax_minor_coverage.set_ylim(0, max(coverage_minor_counts_all)*1.1)


#ax_freq.set_xlabel('Minor allele frequency', fontsize=12)
#ax_freq.set_ylabel('Probability density', fontsize=12)

ax_total_coverage.set_xlabel('Total depth of coverage', fontsize=12)
ax_total_coverage.set_ylabel('Probability density', fontsize=12)

ax_minor_coverage.set_xlabel('Minor allele coverage', fontsize=12)
ax_minor_coverage.set_ylabel('Probability density', fontsize=12)

ax_scatter.set_xlabel('Total depth of coverage', fontsize=12)
ax_scatter.set_ylabel('Median minor allele coverage', fontsize=12)


fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# dpi = 600
fig.savefig("%scoverage_dist.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
plt.close()
