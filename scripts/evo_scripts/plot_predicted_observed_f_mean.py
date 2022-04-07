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



#prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()
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

            predicted_f_mean = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['predicted_f_mean_mapgd']
            predicted_f_mean = numpy.asarray(predicted_f_mean)

            observed_f_mean = prevalence_dict_mapgd[row][clade_type][pi_type][variant_type]['f_mean_mapgd']
            observed_f_mean = numpy.asarray(observed_f_mean)

            predicted_f_mean_no_zeros = predicted_f_mean[(observed_f_mean>0) & (predicted_f_mean>0)]
            observed_f_mean_no_zeros = observed_f_mean[(observed_f_mean>0) & (predicted_f_mean>0)]

            if len(predicted_f_mean_no_zeros) == 0:
                continue

            #all_ = numpy.concatenate([predicted_f_mean_no_zeros,observed_f_mean_no_zeros])

            # Calculate the point density
            #xy = numpy.vstack([observed_f_mean_no_zeros, predicted_f_mean_no_zeros])
            #z = gaussian_kde(xy)(xy)

            # Sort the points by density, so that the densest points are plotted last
            #idx = z.argsort()
            #x, y, z = observed_f_mean_no_zeros[idx], predicted_f_mean_no_zeros[idx], z[idx]
            #ax.scatter(x, y, c=z, cmap=prevalence_utils.variant_cmap_dict[variant_type], s=30, alpha=0.5, edgecolor='', zorder=1)

            sorted_plot_data = prevalence_utils.plot_color_by_pt_dens(observed_f_mean_no_zeros, predicted_f_mean_no_zeros, radius=prevalence_utils.color_radius, loglog=1)
            x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]

            if len(sorted_plot_data[:, 0]) > prevalence_utils.n_points:
                idx_ = numpy.random.choice(len(sorted_plot_data[:, 0]), size=prevalence_utils.n_points, replace=False)
                x = x[idx_]
                y = y[idx_]
                z = z[idx_]

            ax.scatter(x, y, c=numpy.sqrt(z), cmap=prevalence_utils.variant_cmap_dict[variant_type], s=90, alpha=0.9, edgecolors='none', zorder=1)

            all_ = numpy.concatenate([x, y])

            max_ = max(all_)*1.1
            min_ = min(all_)*0.8

            ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2, label='1:1')
            ax.set_xlim([min_, max_])
            ax.set_ylim([min_, max_])

            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)
            ax.set_title(figure_utils.get_pretty_species_name(row), fontsize=16, fontweight='bold', color='k' )

            if (row_idx == 0) and (column_idx == 0):
                ax.legend(loc="upper left", fontsize=14)

            if (row_idx == 0):
                ax.set_ylabel('Predicted mean within-host\nallele frequency across hosts', fontsize=13)

            if column_idx == len(nested_species_list)-1:
                ax.set_xlabel('Observed mean within-host\nallele frequency across hosts', fontsize=13)

            if column_idx == len(nested_species_list)-2:
                if row_idx > len(nested_species_list[-1])-1:
                    ax.set_xlabel('Observed mean within-host\nallele frequency across hosts', fontsize=13)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%spredicted_observed_f_mean_mapgd_%s.pdf" % (config.analysis_directory, variant_type), format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



if __name__=='__main__':

    for variant_type in ['4D', '1D']:

        make_plot(variant_type=variant_type)
