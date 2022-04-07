from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

#import diversity_utils
#import figure_utils
#import parse_midas_data

from scipy.stats import gamma, gaussian_kde

#import calculate_predicted_prevalence_mapgd

from sklearn import svm
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib.lines import Line2D


clade_types = ['all','largest_clade', 'no_strains', 'just_strains']
variant_types = ['4D', '1D']
pi_type = 'pi_include_boundary'


numpy.seterr(divide = 'ignore')



def load_predicted_prevalence_dict_all(test=False):

    intermediate_path = "%spredicted_prevalence_dict_mapgd/" % (config.data_directory)

    prevalence_dict = {}

    for filename in os.listdir(intermediate_path):
        if filename.endswith(".dat"):
            file_path = os.path.join(intermediate_path, filename)
            filename = filename.split('.')[0]

            if filename.count('_') == 3:
                split = filename.rsplit('_', 1)
                species = split[0]
                clade_type = split[1]
            else:
                last_delim = 2
                split = filename.rsplit('_', 2)
                species = split[0]
                clade_type = split[1] + '_' + split[2]

            if test == True:
                if species != 'Bacteroides_ovatus_58035':
                    continue

            if species not in prevalence_dict:
                prevalence_dict[species] = {}

            with open(file_path, 'rb') as handle:
                #b = pickle.load(handle)
                b = pickle.load(handle, encoding='latin1')

            #return b
            prevalence_dict[species][clade_type] = b

    return prevalence_dict




def run_svm():

    f_mean_truncated_all, beta_truncated_all, error_ratio, f_mean_true_all, beta_true_all, f_mean_abc_all, beta_abc_all, predicted_prevalence_truncated_all, observed_prevalence_truncated_all, predicted_prevalence_abc_all, observed_prevalence_abc_all = get_simulated_data()

    log_error_ratio = numpy.log10(error_ratio)

    # run SVM
    y = numpy.copy(log_error_ratio)
    #y[y<=0] = 0
    #y[y>0] = 1

    y[y>=0] = False
    y[y<0] = True

    # 0 = dont run ABC
    # 1 = run ABC

    f_mean_truncated_all_log10 = numpy.log10(f_mean_truncated_all)
    beta_truncated_all_log10 = numpy.log10(beta_truncated_all)

    X = numpy.stack([f_mean_truncated_all_log10, beta_truncated_all_log10], axis=1)

    # select random observations
    rndm_idx = numpy.random.choice(len(y), size=10000, replace=False)

    X_rndm = X[rndm_idx,:]
    y_rndm = y[rndm_idx]

    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

    clf = svm.SVC(kernel='rbf', gamma='auto')
    clf.fit(X_rndm, y_rndm)

    return clf





def get_simulated_data():

    directory = "%sgamma_simulations" % config.data_directory

    f_mean_truncated_all = []
    beta_truncated_all = []

    predicted_prevalence_truncated_all = []
    observed_prevalence_truncated_all = []

    predicted_prevalence_abc_all = []
    observed_prevalence_abc_all = []

    f_mean_true_all = []
    beta_true_all = []
    f_mean_abc_all = []
    beta_abc_all = []

    for filename in os.listdir(directory):

        if filename.endswith(".txt") == False:
            continue

        #0: f_mean
        #1: beta
        #2: predicted_prevalence
        #3: observed_prevalence
        #4: f_mean_truncated
        #5: beta_truncated
        #6: predicted_prevalence_truncated
        #7: observed_prevalence_truncated
        #8: f_mean_abc
        #9: beta_abc
        #10: predicted_prevalence_abc
        #11: observed_prevalence_abcs

        f_mean_true = float(filename.split('_')[2])
        beta_true = float(filename.split('_')[3])

        filepath = os.path.join(directory, filename)
        file = open(filepath, 'r')
        header = file.readline()
        for line in file:

            line = line.strip().split(',')

            line = [float(l) for l in line]

            f_mean_truncated_all.append(line[4])
            beta_truncated_all.append(line[5])

            predicted_prevalence_truncated_all.append(line[6])
            observed_prevalence_truncated_all.append(line[7])

            predicted_prevalence_abc_all.append(line[10])
            observed_prevalence_abc_all.append(line[11])

            f_mean_true_all.append(line[0])
            beta_true_all.append(line[1])

            f_mean_abc_all.append(line[8])
            beta_abc_all.append(line[9])

        file.close()

        filename_split = filename.split('_')
        f_mean = float(filename_split[2])
        beta = float(filename_split[3])
        D = int(filename_split[4].split('.')[0])



        #print("new")
        #print(f_mean_true, )
        #print(beta_true, numpy.mean(beta_abc_file-beta_true)/beta_true)

        #f_mean_true_all.append(f_mean_true)
        #beta_true_all.append(beta_true)
        #f_mean_true_bias_all.append(numpy.mean(f_mean_abc_file-f_mean_true)/f_mean_true)
        #beta_true_bias_all.append(numpy.mean(f_mean_abc_file-f_mean_true)/f_mean_true)


    f_mean_truncated_all = numpy.asarray(f_mean_truncated_all)
    beta_truncated_all = numpy.asarray(beta_truncated_all)
    predicted_prevalence_truncated_all = numpy.asarray(predicted_prevalence_truncated_all)
    observed_prevalence_truncated_all = numpy.asarray(observed_prevalence_truncated_all)
    predicted_prevalence_abc_all = numpy.asarray(predicted_prevalence_abc_all)
    observed_prevalence_abc_all = numpy.asarray(observed_prevalence_abc_all)

    f_mean_true_all = numpy.asarray(f_mean_true_all)
    beta_true_all = numpy.asarray(beta_true_all)
    f_mean_abc_all = numpy.asarray(f_mean_abc_all)
    beta_abc_all = numpy.asarray(beta_abc_all)

    idx = (predicted_prevalence_truncated_all>0) & (observed_prevalence_truncated_all>0) & (predicted_prevalence_abc_all>0) & (observed_prevalence_abc_all>0) & (f_mean_truncated_all>0) & (beta_truncated_all>0)

    f_mean_truncated_all = f_mean_truncated_all[idx]
    beta_truncated_all = beta_truncated_all[idx]
    predicted_prevalence_truncated_all = predicted_prevalence_truncated_all[idx]
    observed_prevalence_truncated_all = observed_prevalence_truncated_all[idx]
    predicted_prevalence_abc_all = predicted_prevalence_abc_all[idx]
    observed_prevalence_abc_all = observed_prevalence_abc_all[idx]

    f_mean_true_all = f_mean_true_all[idx]
    beta_true_all = beta_true_all[idx]
    f_mean_abc_all = f_mean_abc_all[idx]
    beta_abc_all = beta_abc_all[idx]

    error_truncated = numpy.absolute(predicted_prevalence_truncated_all - observed_prevalence_truncated_all)/observed_prevalence_truncated_all
    error_abc = numpy.absolute(predicted_prevalence_abc_all - observed_prevalence_abc_all)/observed_prevalence_abc_all


    #log_error_ratio = numpy.log10(error_abc/error_truncated)
    error_ratio = error_abc/error_truncated
    # remove nans
    idx_error = (numpy.isfinite(error_ratio)) & (beta_truncated_all<=1)


    f_mean_truncated_all = f_mean_truncated_all[idx_error]
    beta_truncated_all = beta_truncated_all[idx_error]
    error_ratio = error_ratio[idx_error]

    # keep simulated values that are in the range of empirical values
    #idx_beta = beta_truncated_all<=1

    #f_mean_truncated_all = f_mean_truncated_all[idx_beta]
    #beta_truncated_all = beta_truncated_all[idx_beta]
    #error_ratio = error_ratio[idx_beta]


    f_mean_true_all = f_mean_true_all[idx_error]
    beta_true_all = beta_true_all[idx_error]
    f_mean_abc_all = f_mean_abc_all[idx_error]
    beta_abc_all = beta_abc_all[idx_error]


    predicted_prevalence_abc_all = predicted_prevalence_abc_all[idx_error]
    observed_prevalence_abc_all = observed_prevalence_abc_all[idx_error]

    predicted_prevalence_truncated_all = predicted_prevalence_truncated_all[idx_error]
    observed_prevalence_truncated_all = observed_prevalence_truncated_all[idx_error]



    return f_mean_truncated_all, beta_truncated_all, error_ratio, f_mean_true_all, beta_true_all, f_mean_abc_all, beta_abc_all, predicted_prevalence_truncated_all, observed_prevalence_truncated_all, predicted_prevalence_abc_all, observed_prevalence_abc_all



def predict_abc_status():

    # classify each site using SVM
    # load the dictionaries to determine whether to run ABC on a given site

    prevalence_dict_mapgd = load_predicted_prevalence_dict_all()
    species_list = list(prevalence_dict_mapgd.keys())

    for species_name in species_list:

        abc_status_dict = {}

        for clade_type in clade_types:

            abc_status_dict[clade_type] = {}

            for variant_type in variant_types:

                f_mean_slm = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['f_mean_slm']
                f_mean_slm = numpy.asarray(f_mean_slm)

                beta_slm = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['beta_slm']
                beta_slm = numpy.asarray(beta_slm)

                sites = prevalence_dict_mapgd[species_name][clade_type][pi_type][variant_type]['sites']

                f_mean_slm_log10 = numpy.log10(f_mean_slm)
                beta_slm_log10 = numpy.log10(beta_slm)

                if len(f_mean_slm_log10) == 0:
                    #print(species_name, clade_type, variant_type)
                    continue

                clf = run_svm()

                X = numpy.asarray([f_mean_slm_log10, beta_slm_log10])
                X = X.T
                yfit = clf.predict(X)

                yfit = numpy.array(yfit, dtype=bool)

                #abc_status_dict[variant_type] = {}

                abc_status_dict[clade_type][variant_type] = {}

                for site_idx, site in enumerate(sites):

                    abc_status_dict[clade_type][variant_type][site] = yfit[site_idx]

                if (clade_type == 'all') and (variant_type == '4D'):

                    print(species_name, sum(yfit)/len(yfit))

        intermediate_filename_template = config.data_directory+"abc_svm_status_dict/%s.dat"
        intermediate_filename = intermediate_filename_template % species_name

        with open(intermediate_filename, 'wb') as handle:
            pickle.dump(abc_status_dict, handle, protocol=2)

        #with open(filepath, mode="wb") as fileObj:
        #    pickle.dump(obj, fileObj, protocol=2)




def plot_svm():

    #f_mean_truncated_all, beta_truncated_all, error_ratio = get_simulated_data()

    f_mean_truncated_all, beta_truncated_all, error_ratio, f_mean_true_all, beta_true_all, f_mean_abc_all, beta_abc_all, predicted_prevalence_truncated_all, observed_prevalence_truncated_all, predicted_prevalence_abc_all, observed_prevalence_abc_all = get_simulated_data()

    log_error_ratio = numpy.log10(error_ratio)


    fig, ax = plt.subplots(figsize=(5,4))

    s = ax.scatter(f_mean_truncated_all, beta_truncated_all, edgecolors='none', s=numpy.absolute(log_error_ratio), alpha=0.6, c=error_ratio, norm=colors.LogNorm(), cmap='RdBu_r')

    cbar = plt.colorbar(mappable=s, ax=ax)
    cbar.set_label('Ratio of prevalence prediction error, ' +  r'$\epsilon_{\mathrm{ABC}}/\epsilon_{\mathrm{Trunc}}$')


    #phase_idx = (numpy.absolute(numpy.log10(error_ratio)) < 0.003)

    #f_mean_truncated_phase = f_mean_truncated_all[phase_idx]
    #beta_truncated_all_phase = beta_truncated_all[phase_idx]

    #ax.scatter(f_mean_truncated_phase, beta_truncated_all_phase, c='k', s=15)


    #print(sum(transition) / len(transition))
    #plt.colorbar()
    #cbar.set_label('Color Intensity')

    ax.set_xlabel('Simulated truncated mean frequency, ' + r'$\bar{f}$', fontsize=11)
    ax.set_ylabel('Simulated truncated squared inverse CV, ' + r'$\beta$', fontsize=11)

    ax.set_xlim([min(f_mean_truncated_all), max(f_mean_truncated_all)])
    ax.set_ylim([min(beta_truncated_all), max(beta_truncated_all)])

    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

    f_mean_truncated_all_log10 = numpy.log10(f_mean_truncated_all)
    beta_truncated_all_log10 = numpy.log10(beta_truncated_all)

    clf = run_svm()

    XX, YY = numpy.mgrid[min(f_mean_truncated_all_log10):max(f_mean_truncated_all_log10):200j, min(beta_truncated_all_log10):max(beta_truncated_all_log10):200j]
    Z = clf.decision_function(numpy.c_[XX.ravel(), YY.ravel()])
    # Put the result into a color plot
    Z = Z.reshape(XX.shape)

    ax.contour(10**XX, 10**YY, Z, colors=["k", "k", "k"], linestyles=["--", "-", "--"], levels=[-0.5, 0, 0.5])

    # '--': support vector
    # '-': separating hyperplane


    legend_elements = [Line2D([0], [0], color='k', lw=2, ls='-', label='Separating hyperplane'),
                       Line2D([0], [0], color='k', lw=2, ls='--', label='Support vector')]

    ax.legend(loc="upper left", handles=legend_elements, fontsize=8)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%sgamma_estimator.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
    plt.close()


plot_svm()

#predict_abc_status()
