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


intermediate_filename_template = config.data_directory+"gamma_simulations/%s.txt"


def run_gamma_simulation(f_mean, beta, coverage, n_hosts, n_iter=1000, min_minor_allele_cov=10, poisson=False):

    # n_iter = number of draws of SNV fluctuation distribution for a given f_mean and beta
    # abc_iter = number of iterations for ABC inference

    sys.stderr.write("Starting simulation...\n")

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

    while n_iter_successes < n_iter:

        f_i_too_many = gamma.rvs(beta, scale=f_mean/beta, size=n_hosts*100)
        f_i_too_many = f_i_too_many[f_i_too_many<1]
        f_i = f_i_too_many[:n_hosts]
        if poisson == True:
            D_i = numpy.random.poisson(lam=coverage, size=n_hosts)
        else:
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

        f_mean_abc_i, beta_abc_i, min_distance = calculate_predicted_prevalence_mapgd.joint_gamma_parameter_simulation(f_mean_truncated_i, beta_truncated_i, D_i)

        # predict prevalence from truncated read counts using mean and beta estimated from *full* distribution
        predicted_prevalence, observed_prevalence, f_mean_, beta_ = calculate_predicted_prevalence_mapgd.predict_prevalence_slm(f_i_truncated, D_i, f_mean=f_mean_i, beta=beta_i, min_minor_allele_cov=10)
        # predict prevalence from truncated read counts using mean and beta estimated from truncated distribution
        predicted_prevalence_truncated, observed_prevalence_truncated, f_mean_truncated_, beta_truncated_ = calculate_predicted_prevalence_mapgd.predict_prevalence_slm(f_i_truncated, D_i, f_mean=f_mean_truncated_i, beta=beta_truncated_i, min_minor_allele_cov=10)
        # predict prevalence from truncated read counts using mean and beta INFERRED using ABC
        predicted_prevalence_abc, observed_prevalence_abc, f_mean_abc, beta_abc = calculate_predicted_prevalence_mapgd.predict_prevalence_slm(f_i_truncated, D_i, f_mean=f_mean_abc_i, beta=beta_abc_i, min_minor_allele_cov=10)

        if any(v == 0 for v in [predicted_prevalence, observed_prevalence, predicted_prevalence_truncated, observed_prevalence_truncated, predicted_prevalence_abc, observed_prevalence_abc]) == True:
            continue

        f_mean_all.append(f_mean_i)
        beta_all.append(beta_i)

        f_mean_truncated_all.append(f_mean_truncated_i)
        beta_truncated_all.append(beta_truncated_i)

        f_mean_abc_all.append(f_mean_abc_i)
        beta_abc_all.append(beta_abc_i)

        predicted_prevalence_all.append(predicted_prevalence)
        observed_prevalence_all.append(observed_prevalence)

        predicted_prevalence_truncated_all.append(predicted_prevalence_truncated)
        observed_prevalence_truncated_all.append(observed_prevalence_truncated)

        predicted_prevalence_abc_all.append(predicted_prevalence_abc)
        observed_prevalence_abc_all.append(observed_prevalence_abc)

        n_iter_successes += 1


    sys.stderr.write("Simulation done!\n")


    f_mean_all = numpy.asarray(f_mean_all)
    beta_all = numpy.asarray(beta_all)

    f_mean_truncated_all = numpy.asarray(f_mean_truncated_all)
    beta_truncated_all = numpy.asarray(beta_truncated_all)

    f_mean_abc_all = numpy.asarray(f_mean_abc_all)
    beta_abc_all = numpy.asarray(beta_abc_all)

    predicted_prevalence_all = numpy.asarray(predicted_prevalence_all)
    observed_prevalence_all = numpy.asarray(observed_prevalence_all)

    predicted_prevalence_truncated_all = numpy.asarray(predicted_prevalence_truncated_all)
    observed_prevalence_truncated_all = numpy.asarray(observed_prevalence_truncated_all)

    predicted_prevalence_abc_all = numpy.asarray(predicted_prevalence_abc_all)
    observed_prevalence_abc_all = numpy.asarray(observed_prevalence_abc_all)

    # error of prevalence from truncated read counts using mean and beta estimated from *full* distribution
    #error_prevalence = numpy.absolute(observed_prevalence_all - predicted_prevalence_all)/observed_prevalence_all

    # error of prevalence from truncated read counts using mean and beta estimated from *truncated* distribution
    #error_prevalence_truncated = numpy.absolute(observed_prevalence_truncated_all - predicted_prevalence_truncated_all)/observed_prevalence_truncated_all

    # error between paramters from full distribution and parameters inferred via ABC
    #error_f_mean_abs = numpy.absolute(f_mean_abc_all - f_mean_all)/f_mean_all
    #error_beta_abs = numpy.absolute(beta_abc_all - beta_all)/beta_all
    # error of prevalence from truncated read counts using mean and beta *inferred* via ABC
    #error_prevalence_abs = numpy.absolute(observed_prevalence_abc_all - predicted_prevalence_abc_all)/observed_prevalence_abc_all
    sys.stderr.write("Saving output......\n")
    record_strs = ["f_mean,beta,predicted_prevalence,observed_prevalence,f_mean_truncated,beta_truncated,predicted_prevalence_truncated,observed_prevalence_truncated,f_mean_abc,beta_abc,predicted_prevalence_abc,observed_prevalence_abc"]

    for idx_ in range(len(f_mean_all)):
        # estimated f_mean,
        record_str = "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g" % (f_mean_all[idx_], beta_all[idx_], predicted_prevalence_all[idx_], observed_prevalence_all[idx_], f_mean_truncated_all[idx_], beta_truncated_all[idx_], predicted_prevalence_truncated_all[idx_], observed_prevalence_truncated_all[idx_], f_mean_abc_all[idx_], beta_abc_all[idx_], predicted_prevalence_abc_all[idx_], observed_prevalence_abc_all[idx_])
        record_strs.append(record_str)


    file_name = 'gamma_simulation_%s_%s_%s' % (str(f_mean), str(beta), str(coverage))
    #file_name = 'test'

    # Write to disk!
    intermediate_filename = intermediate_filename_template % file_name
    #output_file = gzip.open(intermediate_filename,"wb")
    output_file = open(intermediate_filename,"w")
    output_file.write("\n".join(record_strs))
    output_file.close()

    sys.stderr.write("Done!\n")




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







if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--f_mean', type=float)
    parser.add_argument('--beta', type=float)
    # default = median D in B. vulgatus
    parser.add_argument('-D', type=int, default=150)
    # default = n_hosts in B. vulgatus
    parser.add_argument('--n_hosts', type=int, default=300)
    parser.add_argument('--n_iter', type=int, default=1000)
    args = parser.parse_args()


    #f_mean_all = numpy.logspace(-3 , -0.3, num=10, base=10.0)
    #beta_all = numpy.logspace(-2.5, 1, num=10, base=10.0)

    #f_mean_all = numpy.logspace(-3 , -0.3, num=20, base=10.0)
    #beta_all = numpy.logspace(-3 , 1, num=20, base=10.0)

    #print(beta_all)

    f_mean=args.f_mean
    beta=args.beta
    D=args.D
    n_hosts=args.n_hosts
    n_iter=args.n_iter

    run_gamma_simulation(f_mean, beta, D, n_hosts, n_iter=n_iter)
