from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import scipy.stats as stats
import scipy.special as special
import scipy.special

import calculate_predicted_prevalence_mapgd


def test_estimator():

    n_hosts = 100
    coverage = 100
    iter = 1000

    #x_mean_true = 0.05
    #beta_true = 1

    x_mean_all = numpy.logspace(-3 , -0.3, num=20, base=10.0)
    beta_all = numpy.logspace(-2.5, 1, num=20, base=10.0)
    #x_mean_all = [0.001]
    #beta_all = [0.0031622776601683794]
    #x_mean_all = 10**x_mean_all_log10
    #beta_all = 10**beta_all_log10

    error_ratio_dict = {}

    for x_mean_true in x_mean_all:

        error_ratio_dict[x_mean_true] = {}

        for beta_true in beta_all:

            print(x_mean_true, beta_true)

            error_prevalence_improvement_all = []
            error_ratio_all = []
            min_distance_slm_all = []
            error_f_mean_all = []
            error_beta_all = []
            count = 0
            #for i in range(iter):
            while count < iter:

                D = numpy.random.poisson(lam=coverage, size=n_hosts)

                x_i = stats.gamma.rvs(beta_true, scale=x_mean_true/beta_true, size=n_hosts)

                x_i_to_keep = x_i[x_i<=1]
                D_i_to_keep = D[x_i<=1]
                A_i = numpy.random.binomial(D_i_to_keep, x_i_to_keep)
                observed_prevalence = sum(A_i>0)/len(A_i)
                A_i[A_i<10] = 0
                observed_prevalence_truncated = sum(A_i>0)/len(A_i)
                f_i = A_i/D_i_to_keep

                observed_x_mean_truncated = numpy.mean(f_i)
                if observed_x_mean_truncated == 0:
                    continue

                observed_beta_truncated = (numpy.mean(f_i)**2)/numpy.var(f_i)

                x_mean_best, beta_best, min_distance_slm = calculate_predicted_prevalence_mapgd.joint_gamma_parameter_simulation(observed_x_mean_truncated, observed_beta_truncated, D_i_to_keep, iter = 10000)

                error_f_mean = numpy.absolute(x_mean_best - observed_x_mean_truncated) / observed_x_mean_truncated
                error_beta = numpy.absolute(beta_best - observed_beta_truncated) / observed_beta_truncated


                predicted_prevalence_slm, observed_prevalence, f_mean, beta = calculate_predicted_prevalence_mapgd.predict_prevalence_slm(f_i, D_i_to_keep, f_mean=observed_x_mean_truncated, beta=observed_beta_truncated, min_minor_allele_cov=10)
                predicted_prevalence_slm_best, observed_prevalence_best, f_mean_best, beta_best = calculate_predicted_prevalence_mapgd.predict_prevalence_slm(f_i, D_i_to_keep, f_mean=x_mean_best, beta=beta_best, min_minor_allele_cov=10)


                error = numpy.absolute(observed_prevalence_truncated - predicted_prevalence_slm)/observed_prevalence_truncated
                error_best = numpy.absolute(observed_prevalence_truncated - predicted_prevalence_slm_best)/observed_prevalence_truncated

                error_prevalence_improvement = (error_best-error)/error
                error_prevalence_improvement_all.append(error_prevalence_improvement)

                error_ratio_all.append(error_best/error)

                min_distance_slm_all.append(min_distance_slm)

                error_f_mean_all.append(error_f_mean)
                error_beta_all.append(error_beta)

                count += 1


            error_ratio_dict[x_mean_true][beta_true] = {}
            error_ratio_dict[x_mean_true][beta_true]['prevalence_improvement'] = numpy.mean(error_prevalence_improvement_all)
            error_ratio_dict[x_mean_true][beta_true]['prevalence_error_ratio'] = numpy.mean(error_ratio_all)
            error_ratio_dict[x_mean_true][beta_true]['min_distance'] = numpy.mean(min_distance_slm_all)
            error_ratio_dict[x_mean_true][beta_true]['error_f_mean'] = numpy.mean(error_f_mean_all)
            error_ratio_dict[x_mean_true][beta_true]['error_beta'] = numpy.mean(error_beta_all)



    path = "%serror_ratio_dict.dat" % config.data_directory

    with open(path, 'wb') as outfile:
        pickle.dump(error_ratio_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)









test_estimator()
