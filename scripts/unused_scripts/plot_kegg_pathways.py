import  matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import parse_midas_data
#import pylab
#from pylab import *
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
import os
#import pandas
import parse_patric
import pickle
import sample_utils
import parse_HMP_data

import itertools

import scipy.stats as stats

#kegg_pi_Eubacterium_rectale_56927.dat


good_species_list = parse_midas_data.parse_good_species_list()
#good_species_list = ["Eubacterium_rectale_56927"]

kegg_dict_all_species = {}

for species_name in good_species_list:


    kegg_path = '%skegg_pi/kegg_pi_%s.dat' % (parse_midas_data.data_directory, species_name)

    if os.path.exists(kegg_path) == False:
        continue

    kegg_dict = pickle.load( open(kegg_path, "rb" ) )

    pathways = kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys()

    for pathway in kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys():

        pathway_array = kegg_dict['fraction_nonsynonymous_per_pathway_core'][pathway]

        n_hosts = pathway_array.shape[0]

        values = pathway_array[numpy.triu_indices(n_hosts, 1)]

        #kegg_dict_all_species[species_name][pathway] = numpy.mean(values[values>0])

        if pathway not in kegg_dict_all_species:
            kegg_dict_all_species[pathway] = []

        if len(values[(values>0) & (values<1)]) <= 20:
            continue

        kegg_dict_all_species[pathway].append(numpy.mean(values[(values>0) & (values<1)]))


#kegg1 = [kegg_dict_all_species[species]['Alanine, aspartate and glutamate metabolism'] for species in kegg_dict_all_species.keys() if 'Alanine, aspartate and glutamate metabolism' in kegg_dict_all_species[species]]
#kegg2 = [kegg_dict_all_species[species]['Lipopolysaccharide biosynthesis'] for species in kegg_dict_all_species.keys() if 'Lipopolysaccharide biosynthesis' in kegg_dict_all_species[species]]

fraction_nonsynonymous_per_pathway_core_means = []

pathway_names = kegg_dict_all_species.keys()

paired_means = [ [pathway_name, numpy.mean(kegg_dict_all_species[pathway_name]), numpy.var(kegg_dict_all_species[pathway_name]) ] for pathway_name in kegg_dict_all_species.keys() if (len(kegg_dict_all_species[pathway_name]) > 10) and (pathway_name != '') ]

paired_means_sorted = sorted(paired_means, key = lambda x: float(x[1]))


#ax_pca = fig.add_subplot(gs[0:3, 2:4])
#fig, ax = plt.subplots(figsize = (12, 9))
fig = plt.figure(figsize = (12, 6))

ax_mean = plt.subplot2grid((1, 2), (0, 0))
ax_taylors = plt.subplot2grid((1, 2), (0, 1))


for paired_means_sorted_i_idx, paired_means_sorted_i in enumerate(paired_means_sorted):

    paired_i_name, paired_i_mean, paired_i_var  = paired_means_sorted_i

    ax_mean.scatter(kegg_dict_all_species[paired_i_name], list(itertools.repeat(paired_means_sorted_i_idx, len(kegg_dict_all_species[paired_i_name]))),
            s = 15, \
            linewidth=1, facecolors='dodgerblue', \
            edgecolors='k', marker='o', \
            alpha=0.5, zorder=1)

    ax_mean.scatter(paired_i_mean, paired_means_sorted_i_idx,
            s = 28, \
            linewidth=1, facecolors='k', \
            edgecolors='k', marker='o', \
            alpha=1, zorder=2)





labels = [x[0] for x in paired_means_sorted]

ax_mean.tick_params(axis='y', which='y',length=0)

ax_mean.set_yticks(numpy.arange(len(labels)))

ax_mean.set_yticklabels(labels)
ax_mean.yaxis.set_tick_params(labelsize=5.5)

ax_mean.set_xlabel('Proportion of nonsynonymous substitutions', fontsize=12)



means = [x[1] for x in paired_means_sorted]
variances = [x[2] for x in paired_means_sorted]


ax_taylors = inset_axes(ax_fmax, width="100%", height="100%", loc='lower right', bbox_to_anchor=(0.12,0.07,0.4,0.38), bbox_transform=ax_fmax.transAxes)

ax_taylors.scatter(means, variances, s = 50, \
        linewidth=1.5, facecolors='dodgerblue', \
        edgecolors='k', marker='o', \
        alpha=0.9)


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(means), numpy.log10(variances))


x_log10_fit_range =  numpy.linspace(numpy.log10(min(means) * 0.5), numpy.log10(max(means) * 1.5), 10000)

y_fit_range = 10 ** (slope*x_log10_fit_range + intercept)
ax_taylors.plot(10**x_log10_fit_range, y_fit_range, c='k', lw=3, linestyle='--', zorder=2)


ax_taylors.text(0.3,0.9, r'$\sigma^{{2}}_{{ P_{{N}} }} \sim \left \langle P_{{N}} \right \rangle^{{{}}}$'.format(str(round(slope, 3)) ), fontsize=11, color='k', ha='center', va='center', transform=ax_taylors.transAxes  )



ax_taylors.set_xscale('log', basex=10)
ax_taylors.set_yscale('log', basey=10)

ax_taylors.set_xlabel('Mean proportion of nonsynonymous substitutions', fontsize=12)
ax_taylors.set_ylabel('Variance of proportion of nonsynonymous substitutions', fontsize=12)



fig.subplots_adjust(hspace=0.4, wspace=0.35) #hspace=0.3, wspace=0.5
fig_name =  "%skegg.png" % parse_midas_data.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
