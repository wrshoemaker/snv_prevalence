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

import parse_midas_data

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.stats import gamma, gaussian_kde

import calculate_predicted_occupancy

data_directory = config.data_directory
allowed_variant_types = set(['1D','4D'])

#color_dict = {'4D': 'b', '1D': 'r'}

#intermediate_filename_template = config.data_directory+"coprevalence_f0/%s.dat"


#good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list = ['Bacteroides_vulgatus_57955']



prevalence_dict = calculate_predicted_occupancy.load_predicted_prevalence_subsample_dict()

good_species_list = prevalence_dict.keys()


def plot_obs_pred():

    for species_name in good_species_list:

        #intermediate_filename_template = config.data_directory+"predicted_observed_prevalence/%s.dat"

        #intermediate_filename = intermediate_filename_template % species_name

        #if os.path.isfile(intermediate_filename) == False:
        #    continue

        #with open(intermediate_filename, 'rb') as handle:
        #    snp_dict = pickle.load(handle)

        #print(species_name)


        predicted_prevalence = snp_dict['4D']['predicted_prevalence']
        predicted_prevalence = numpy.asarray(predicted_prevalence)

        observed_prevalence = snp_dict['4D']['observed_prevalence']
        observed_prevalence = numpy.asarray(observed_prevalence)

        f_max = snp_dict['4D']['f_max']
        f_max = numpy.asarray(f_max)

        f_mean = snp_dict['4D']['f_mean_no_zeros']
        f_var = snp_dict['4D']['f_var_no_zeros']

        f_mean = numpy.asarray(f_mean)
        f_var = numpy.asarray(f_var)

        #f_mean = f_mean[f_var>0]
        #f_var = f_var[f_var>0]



        #fig = plt.figure(figsize = (4,4)) #

        fig, ax = plt.subplots(figsize=(4,4))

        #ax_obs_pred = plt.subplot2grid((1, 4), (0, 0), colspan=1)

        #ax_obs_pred = plt.subplot2grid((1, 4), (0, 0), colspan=1)

        #ax_fmax = plt.subplot2grid((1, 4), (0, 1), colspan=1)

        #ax_fmax_prev = plt.subplot2grid((1, 4), (0, 2), colspan=1)

        #ax_f_mean_var = plt.subplot2grid((1, 4), (0, 3), colspan=1)


        # Calculate the point density
        xy = numpy.vstack([observed_prevalence,predicted_prevalence])
        z = gaussian_kde(xy)(xy)

        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = observed_prevalence[idx], predicted_prevalence[idx], z[idx]

        ax_ = ax.scatter(x, y, c=z, cmap=cm.Blues, s=40, alpha=0.9, edgecolor='')
        plt.colorbar(ax_)


        #ax_obs_pred.scatter(observed_prevalence, predicted_prevalence, alpha=0.1)

        ax.plot([min(observed_prevalence)*0.8, max(observed_prevalence)*1.1],[min(observed_prevalence)*0.8, max(observed_prevalence)*1.1], ls='--', c='k')

        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)

        ax.set_xlabel('Observed SNP prevalence', fontsize=12)
        ax.set_ylabel('Predicted SNP prevalence, gamma', fontsize=12)

        #ax.set_title(species_name)


        ax.set_title(species_name)


        #print(numpy.mean(numpy.abs(predicted_prevalence-observed_prevalence)))

        #ax_fmax.scatter(f_mean, predicted_prevalence-observed_prevalence, alpha=0.1)

        #ax_fmax.set_xlabel('Mean frequency', fontsize=12)
        #ax_fmax.set_ylabel('Predicted - observed prevalence', fontsize=12)


        #ax_fmax.set_xscale('log', basex=10)

        #ax_fmax.axhline(y=0, c='k', ls='--')

        #ax_fmax.set_xlim([0.0001, 1])



        #ax_fmax_prev.scatter(f_mean, observed_prevalence, alpha=0.1)
        #ax_fmax_prev.set_xlabel('Mean frequency', fontsize=12)
        #ax_fmax_prev.set_ylabel('Observed prevalence', fontsize=12)
        #ax_fmax_prev.set_xscale('log', basex=10)
        #ax_fmax_prev.set_yscale('log', basey=10)

        #ax_fmax_prev.set_xlim([0.0001, 1])
        #ax_fmax_prev.set_ylim([0.01, 1])

        #ax_f_mean_var.scatter(f_mean, f_var, alpha = 0.1)
        #ax_f_mean_var.set_xscale('log', basex=10)
        #ax_f_mean_var.set_yscale('log', basey=10)

        #ax_f_mean_var.set_xlim([0.00001, 0.6])
        #ax_f_mean_var.set_ylim([0.000001, 0.6])


        fig.tight_layout()
        fig.savefig("%spredict_occupancy_snps/%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()



def plot_obs_pred_strain_comparison():

    for species_name in good_species_list:

        print(species_name)

        fig = plt.figure(figsize = (8,4))

        directories_label_dict = {'predicted_observed_prevalence': 'All samples', 'predicted_observed_prevalence_no_strains': 'Samples without strain structure'}

        directories = ['predicted_observed_prevalence', 'predicted_observed_prevalence_no_strains']
        #directories = ['predicted_observed_prevalence']

        for directory_idx, directory in enumerate(directories):

            #intermediate_filename_template = config.data_directory+directory+"/%s.dat"
            #intermediate_filename = intermediate_filename_template % species_name
            #if os.path.isfile(intermediate_filename) == False:
            #    continue

            #with open(intermediate_filename, 'rb') as handle:
            #    snp_dict = pickle.load(handle)

            if directory not in prevalence_dict[species_name]:
                continue


            predicted_prevalence = prevalence_dict[species_name][directory]['predicted_prevalence_to_plot']
            predicted_prevalence = numpy.asarray(predicted_prevalence)



            observed_prevalence = prevalence_dict[species_name][directory]['observed_prevalence_to_plot']
            observed_prevalence = numpy.asarray(observed_prevalence)


            mae = prevalence_dict[species_name][directory]['MAE']


            #subset_idx = numpy.random.randint(0, high=len(predicted_prevalence), size=10000)
            subset_idx = numpy.random.choice(len(predicted_prevalence), size=10000, replace=False)



            predicted_prevalence_subset = predicted_prevalence[subset_idx]
            observed_prevalence_subset = observed_prevalence[subset_idx]


            ax = plt.subplot2grid((1, 2), (0, directory_idx), colspan=1)

            # Calculate the point density
            xy = numpy.vstack([observed_prevalence_subset,predicted_prevalence_subset])
            z = gaussian_kde(xy)(xy)

            # Sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            x, y, z = observed_prevalence_subset[idx], predicted_prevalence_subset[idx], z[idx]


            ax_ = ax.scatter(x, y, c=z, cmap=cm.Blues, s=40, alpha=0.9, edgecolor='')
            plt.colorbar(ax_)



            ax.plot([min(observed_prevalence)*0.8, max(observed_prevalence)*1.1],[min(observed_prevalence)*0.8, max(observed_prevalence)*1.1], ls='--', c='k')

            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)

            ax.set_xlabel('Observed SNP prevalence', fontsize=12)
            ax.set_ylabel('Predicted SNP prevalence, gamma', fontsize=12)

            #ax.set_title(species_name)

            ax.set_title(directories_label_dict[directory])


        fig.tight_layout()
        fig.savefig("%spredict_occupancy_snps_comparison/%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()


plot_obs_pred_strain_comparison()
