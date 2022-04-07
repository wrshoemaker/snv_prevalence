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

import calculate_predicted_prevalence_mapgd


numpy.random.seed(123456789)


coverage = 150
n_hosts =100
f_mean = 0.1
beta = 0.2
min_minor_allele_cov=10

n_iter = 50

coverage_null = numpy.asarray([coverage]*n_hosts)

n_iter_successes = 0

#for i in range(n_iter):

f_mean_all = []
beta_all = []

f_mean_truncated_all = []
beta_truncated_all = []

f_mean_abc_all = []
beta_abc_all = []

predicted_prevalence_all = []
observed_prevalence_all = []

predicted_prevalence_truncated_all = []
observed_prevalence_truncated_all = []

predicted_prevalence_abc_all = []
observed_prevalence_abc_all = []

error_euc_all = []
error_man_all = []

while n_iter_successes < n_iter:

    f_i_too_many = gamma.rvs(beta, scale=f_mean/beta, size=n_hosts*100)
    f_i_too_many = f_i_too_many[f_i_too_many<1]
    f_i = f_i_too_many[:n_hosts]
    D_i = coverage_null

    A_i = numpy.random.binomial(D_i, f_i)
    f_i_naive = A_i/D_i

    f_mean_i = numpy.mean(f_i_naive)
    f_var_i = numpy.var(f_i_naive)

    if (f_mean_i == 0) or (f_var_i == 0):
        continue

    # truncate read counts
    A_i[A_i<min_minor_allele_cov] = 0
    f_i_truncated = A_i/D_i
    f_mean_truncated_i = numpy.mean(f_i_truncated)
    f_var_truncated_i = numpy.var(f_i_truncated)

    if (f_mean_truncated_i == 0) or (f_var_truncated_i == 0):
        continue

    beta_i = (f_mean_i**2) / f_var_i
    beta_truncated_i = (f_mean_truncated_i**2) / f_var_truncated_i

    f_mean_abc_i, beta_abc_i, min_distance = calculate_predicted_prevalence_mapgd.joint_gamma_parameter_simulation(f_mean_truncated_i, beta_truncated_i, D_i, distance='euclidean', iter = 10000)
    f_mean_abc_i_manhattan, beta_abc_i_manhattan, min_distance_manhattan = calculate_predicted_prevalence_mapgd.joint_gamma_parameter_simulation(f_mean_truncated_i, beta_truncated_i, D_i, distance='manhattan', iter = 10000)

    error_euc = numpy.absolute(beta_i - beta_abc_i)/beta_i
    error_man = numpy.absolute(beta_i - beta_abc_i_manhattan)/beta_i


    error_euc_all.append(error_euc)
    error_man_all.append(error_man)

    n_iter_successes += 1


print(numpy.mean(error_euc_all), numpy.mean(error_man_all) )
