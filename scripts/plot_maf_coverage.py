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

from scipy.stats import gamma
import scipy.special
from scipy.integrate import quad


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

import calculate_predicted_prevalence_mapgd


species_to_run = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                    'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                    'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035',
                    'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318',
                    'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208',
                    'Dialister_invisus_61905', 'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985',
                    'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']



max_cov = 20
range_ = list(range(max_cov+1))


fig, ax = plt.subplots(figsize=(4.5,4))
#gs = gridspec.GridSpec(nrows=1, ncols=2)
#fig = plt.figure(figsize = (8, 4))

#ax_freq = fig.add_subplot(gs[0, 0])
#ax_coverage = fig.add_subplot(gs[0, 1])

for species_name in species_to_run:


    intermediate_filename_template = config.data_directory+"mapgd_output_dicts/%s.dat"
    intermediate_filename = intermediate_filename_template % species_name

    if os.path.isfile(intermediate_filename) == False:
        continue

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)

    frequency_minor_allele = b['frequency_minor_allele']
    coverage_minor_allele = b['coverage_minor_allele']

    coverage_minor_allele = [int(round(c)) for c in coverage_minor_allele]

    frequency_minor_allele = numpy.asarray(frequency_minor_allele)
    coverage_minor_allele = numpy.asarray(coverage_minor_allele)

    frequency_minor_allele = numpy.log10(frequency_minor_allele)

    #ax_freq.hist(frequency_minor_allele, bins=30, density=True, histtype='step')

    coverage_minor_allele_counts = numpy.asarray([sum(coverage_minor_allele==i) for i in range(max_cov+1)])
    coverage_minor_allele_counts = coverage_minor_allele_counts/len(coverage_minor_allele)

    #ax_coverage.hist(coverage_minor_allele, bins=30, density=True, histtype='step')

    ax.scatter(range_, coverage_minor_allele_counts)



ax.set_xlim(0, 20)

#ax.set_ylim(0, 0.18)

#ax_freq.set_xlabel('Minor allele frequency', fontsize=12)
#ax_freq.set_ylabel('Probability density', fontsize=12)

ax.set_xlabel('MAPGD minor allele read count', fontsize=12)
ax.set_ylabel('Probability density', fontsize=12)

fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
# dpi = 600
fig.savefig("%sminor_allele_dist.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
plt.close()
