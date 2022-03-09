from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_prevalence_mapgd

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

#clade_types = ['all','largest_clade', 'no_strains']


prevalence_dict_mapgd = calculate_predicted_prevalence_mapgd.load_predicted_prevalence_dict_all()
species_list = list(prevalence_dict_mapgd.keys())
species_list.sort()
n_species_row=5
nested_species_list = [species_list[x:x+n_species_row] for x in range(0, len(species_list), n_species_row)]

clade_type = 'all'
pi_type = 'pi_include_boundary'


def make_plot(variant_type):

    gs = gridspec.GridSpec(nrows=len(nested_species_list), ncols=2*n_species_row)
    fig = plt.figure(figsize = (40, 20))

    for column_idx, column in enumerate(nested_species_list):

        for row_idx, row in enumerate(column):

            ax = fig.add_subplot(gs[column_idx, row_idx])

            predicted_prevalence_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_prevalence_mapgd']
            predicted_prevalence_slm = numpy.asarray(predicted_prevalence_slm)

            observed_prevalence_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['observed_prevalence_mapgd']
            observed_prevalence_slm = numpy.asarray(observed_prevalence_slm)

            f_mean_slm = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['f_max_mapgd']
            f_mean_slm = numpy.asarray(f_mean_slm)

            sites = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['sites']

            idx_ = (observed_prevalence_slm>0) & (predicted_prevalence_slm>0) & (f_mean_slm>0)

            predicted_prevalence_slm = predicted_prevalence_slm[idx_]
            observed_prevalence_slm = observed_prevalence_slm[idx_]
            f_mean_slm = f_mean_slm[idx_]

            f_mean_slm_log10 = numpy.log10(f_mean_slm)
            # relationship between f_mean and prevalence
            hist, bin_edges = numpy.histogram(f_mean_slm_log10, density=True, bins=40)
            bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]
            prevalence_predicted_mean = []
            for i in range(0, len(bin_edges)-1):
                idx_i = (f_mean_slm_log10 > bin_edges[i]) & (f_mean_slm_log10 <= bin_edges[i+1])
                prevalence_predicted_mean.append(numpy.mean(numpy.log10(predicted_prevalence_slm[idx_i])))
            bins_mean = numpy.asarray(bins_mean)
            prevalence_predicted_mean = numpy.asarray(prevalence_predicted_mean)

            #all_ = numpy.concatenate([f_mean_slm, observed_prevalence_slm])
            #xy = numpy.vstack([f_mean_slm, observed_prevalence_slm])
            #z = gaussian_kde(xy)(xy)
            # Sort the points by density, so that the densest points are plotted last
            #idx = z.argsort()
            #x, y, z = f_mean_slm[idx], observed_prevalence_slm[idx], z[idx]
            #ax.scatter(x, y, c=z, cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=1, edgecolor='', zorder=1, label = 'Observed')

            sorted_plot_data = prevalence_utils.plot_color_by_pt_dens(f_mean_slm, observed_prevalence_slm, radius=prevalence_utils.color_radius, loglog=1)
            x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]

            if len(sorted_plot_data[:, 0]) > prevalence_utils.n_points:
                idx_ = numpy.random.choice(len(sorted_plot_data[:, 0]), size=prevalence_utils.n_points, replace=False)
                x = x[idx_]
                y = y[idx_]
                z = z[idx_]

            ax.scatter(x, y, c=numpy.sqrt(z), cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolors='none', zorder=1)
            all_ = numpy.concatenate([x, y])

            ax.plot(10**bins_mean, 10**prevalence_predicted_mean, c='k', lw = 4, ls='--', label = 'Predicted', zorder=2)

            ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=12, fontweight='bold', color='k' )
            ax.set_xlim([min(f_mean_slm), max(f_mean_slm)])
            #ax.set_ylim([min(observed_prevalence_slm), max(observed_prevalence_slm)])
            ax.set_ylim([min(predicted_prevalence_slm), max(observed_prevalence_slm)])

            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)

            if (row_idx == 0) and (column_idx == 0):
                ax.legend(loc="upper left", fontsize=14)

            if (row_idx == 0):
                ax.set_ylabel('Allele prevalence', fontsize=14)


            if column_idx == len(nested_species_list)-1:
                ax.set_xlabel('Max. within-host allele\nfrequency across hosts, ' + r'$f_{max}$', fontsize=14)


            if column_idx == len(nested_species_list)-2:
                if row_idx > len(nested_species_list[-1])-1:
                    ax.set_xlabel('Max. within-host allele\nfrequency across hosts, ' + r'$f_{max}$', fontsize=14)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%sf_max_vs_prevalence_%s.pdf" % (config.analysis_directory, variant_type), format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
    plt.close()



if __name__=='__main__':

    for variant_type in ['4D', '1D']:

        make_plot(variant_type)
