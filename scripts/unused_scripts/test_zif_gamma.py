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

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import calculate_predicted_prevalence_mapgd
#import calculate_predicted_prevalence

numpy.random.seed(123456789)

# low to high values
colors = ['lightskyblue', 'dodgerblue', 'darkblue']

beta = 0.8

def simulate_zif_gamma(beta, x_mean, prob_absence, cov=50, n_hosts=100):

    shape = beta
    rate = beta/x_mean
    scale = 1/rate

    n_hosts = 100

    # 50 fold coverage
    #D = numpy.repeat(50, n_hosts)
    D = numpy.random.poisson(lam=cov, size=n_hosts)

    n_absences = numpy.random.binomial(n_hosts, prob_absence)

    x = stats.gamma.rvs(shape, scale=scale, size=n_hosts-n_absences)

    x_zif = numpy.concatenate((x, numpy.asarray([0]*n_absences)))

    x_zif_to_keep = x_zif[x_zif<=1]
    D_to_keep = D[x_zif<=1]

    A = numpy.random.binomial(D_to_keep, x_zif_to_keep)

    observed_prevalence = sum(A>0)/len(A)

    predicted_prevalence = 1 - numpy.mean((1 + (x_mean/beta)*D_to_keep)**(-1*beta))

    if observed_prevalence > 0:

        #relative_error = numpy.absolute(predicted_prevalence - observed_prevalence)/observed_prevalence
        relative_error = (predicted_prevalence - observed_prevalence)/observed_prevalence

    else:
        relative_error = float("nan")

    return observed_prevalence, predicted_prevalence, relative_error




def error_plot():


    #x_mean = 0.001
    prob_absence_all = numpy.linspace(0, 0.9, num=20)

    c_all = ['k', 'r', 'b']

    fig, ax = plt.subplots(figsize=(4,4))

    for x_mean_i_idx, x_mean_i in enumerate([0.001, 0.01, 0.1]):

        mean_relative_error_all = []

        print(x_mean_i)

        for prob_absence_i in prob_absence_all:

            relative_error_all = []

            for i in range(1000):

                observed_prevalence, predicted_prevalence, relative_error = simulate_zif_gamma(beta, x_mean_i, prob_absence_i)
                relative_error_all.append(relative_error)

            relative_error_all = numpy.asarray(relative_error_all)

            mean_relative_error_all.append(numpy.mean(relative_error_all[~numpy.isnan(relative_error_all)]))

        mean_relative_error_all = numpy.asarray(mean_relative_error_all)

        ax.scatter(prob_absence_all, mean_relative_error_all, alpha=1, c=c_all[x_mean_i_idx], label='mean SNV frequency = %0.3f' % x_mean_i )#, c='#87CEEB')




    ax.set_title("Error in prevalence predictions when\nthe true distribution is a zero-inflated gamma", fontsize=12, fontweight='bold' )

    ax.set_xlabel('Probability that a SNV is absent', fontsize=12)
    ax.set_ylabel('(Observed - predicted)/observed', fontsize=12)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%stest_zif_gamma.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()




def plot_observed_vs_predicted_zif_gamma():

    x_mean_all = numpy.random.uniform(low=-3, high=-1, size=1000)
    x_mean_all = 10**x_mean_all

    observed_prevalence_all = []
    predicted_prevalence_all = []
    for x_mean_i in x_mean_all:

        observed_prevalence, predicted_prevalence, relative_error = simulate_zif_gamma(beta, x_mean_i, 0.5)

        if (observed_prevalence>0) and (predicted_prevalence>0):
            observed_prevalence_all.append(observed_prevalence)
            predicted_prevalence_all.append(predicted_prevalence)

    observed_prevalence_all = numpy.asarray(observed_prevalence_all)
    predicted_prevalence_all = numpy.asarray(predicted_prevalence_all)


    fig, ax = plt.subplots(figsize=(4,4))

    ax.scatter(observed_prevalence_all, predicted_prevalence_all, alpha=0.7)

    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)


    all_ = numpy.concatenate([observed_prevalence_all,predicted_prevalence_all])

    max_ = max(all_)*1.1
    min_ = min(all_)*0.8

    #print(min(f_max_no_zeros))

    ax.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
    ax.set_xlim([min_, max_])
    ax.set_ylim([min_, max_])

    ax.set_xlabel('Observed prevalence', fontsize=12)
    ax.set_ylabel('Predicted prevalence', fontsize=12)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%stest_predicted_prevalence_zif_gamma.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()





def test_cumulative_prob_absence():

    min_minor_allele_cov = 10
    f_mean = 0.01
    f_beta = 1

    depths = numpy.asarray([200]*10)

    prob_zero = (1+ ((f_mean/f_beta)*depths))**(-1*f_beta )
    predicted_prevalence_slm = 1 - numpy.mean(prob_zero)
    #print(predicted_prevalence_slm)

    if min_minor_allele_cov > 1:
        # range() excludes upper bound
        cumulative_sum = 0
        for n_i in range(0, min_minor_allele_cov):
            cumulative_sum += (scipy.special.gamma(f_beta+n_i) / numpy.math.factorial(n_i)) * (((f_mean*depths)/ (f_beta + f_mean*depths)) ** n_i) / scipy.special.gamma(f_beta)

        prob_zero *= cumulative_sum


    predicted_prevalence_slm = 1 - numpy.mean(prob_zero)



def test_gamma_parameter_sampling():

    iter = 1000

    x_mean_all = numpy.logspace(-2, -0.4 , num=100)
    beta_all = numpy.logspace(-2, 1 , num=100)

    n_hosts = 100
    coverage = 500


    fig = plt.figure(figsize = (8, 4)) #
    fig.subplots_adjust(bottom= 0.15)

    ax_mean = plt.subplot2grid((1, 2), (0, 0), colspan=1)
    ax_beta = plt.subplot2grid((1, 2), (0, 1), colspan=1)

    betas = [0.1, 1, 10]
    x_means = [0.01, 0.03, 0.1]

    # inverse transform sampling

    ax_mean_x_all = []
    ax_mean_y_all = []
    for beta_i_idx, beta_i in enumerate(betas):

        x_mean_estimate_mean_all = []
        for x_mean in x_mean_all:

            shape = beta_i
            rate = beta_i/x_mean
            scale = 1/rate

            x_mean_estimate_all = []
            for i in range(iter):

                # 50 fold coverage
                D = numpy.random.poisson(lam=coverage, size=n_hosts)

                x = stats.gamma.rvs(shape, scale=scale, size=n_hosts)
                x_to_keep = x[x<=1]
                D_to_keep = D[x<=1]

                A = numpy.random.binomial(D_to_keep, x_to_keep)
                A[A < 10] = 0
                freqs = A/D_to_keep
                x_mean_estimate = numpy.mean(freqs)

                x_mean_estimate_all.append(x_mean_estimate)

            x_mean_estimate_mean = numpy.mean(x_mean_estimate_all)

            x_mean_estimate_mean_all.append(x_mean_estimate_mean)


        ax_mean.scatter(x_mean_all, x_mean_estimate_mean_all, c=colors[beta_i_idx], label=r'$\beta=$' + str(round(beta_i, 3)), alpha=0.7)

        ax_mean_x_all.extend(x_mean_all)
        ax_mean_y_all.extend(x_mean_estimate_mean_all)



    ax_beta_x_all = []
    ax_beta_y_all = []
    for x_mean_i_idx, x_mean_i in enumerate(x_means):

        beta_estimate_mean_all = []
        for beta in beta_all:

            shape = beta
            rate = beta/x_mean_i
            scale = 1/rate

            beta_estimate_all = []
            for i in range(iter):

                # 50 fold coverage
                D = numpy.random.poisson(lam=coverage, size=n_hosts)

                x = stats.gamma.rvs(shape, scale=scale, size=n_hosts)
                x_to_keep = x[x<=1]
                D_to_keep = D[x<=1]

                A = numpy.random.binomial(D_to_keep, x_to_keep)
                A[A < 10] = 0
                freqs = A/D_to_keep
                beta_estimate = (numpy.mean(freqs)**2) / numpy.var(freqs)

                beta_estimate_all.append(beta_estimate)

            beta_estimate_mean = numpy.mean(beta_estimate_all)

            beta_estimate_mean_all.append(beta_estimate_mean)


        ax_beta.scatter(beta_all, beta_estimate_mean_all, c=colors[x_mean_i_idx], label=r'$\bar{x}=$' + str(round(x_mean_i, 3)), alpha=0.7)

        ax_beta_x_all.extend(beta_all)
        ax_beta_y_all.extend(beta_estimate_mean_all)



    ax_mean.set_xscale('log', basex=10)
    ax_mean.set_yscale('log', basey=10)
    all_ax_mean = numpy.concatenate([ax_mean_x_all, ax_mean_y_all])
    max_ = max(all_ax_mean)*1.1
    min_ = min(all_ax_mean)*0.8
    ax_mean.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
    ax_mean.set_xlim([min_, max_])
    ax_mean.set_ylim([min_, max_])
    ax_mean.set_xlabel('Real mean frequency', fontsize=12)
    ax_mean.set_ylabel('Estimate of mean frequency', fontsize=12)
    ax_mean.legend(loc="upper left", fontsize=8)


    ax_beta.set_xscale('log', basex=10)
    ax_beta.set_yscale('log', basey=10)
    all_ax_beta = numpy.concatenate([ax_beta_x_all, ax_beta_y_all])
    max_ = max(all_ax_beta)*1.1
    min_ = min(all_ax_beta)*0.8
    ax_beta.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
    ax_beta.set_xlim([min_, max_])
    ax_beta.set_ylim([min_, max_])
    ax_beta.set_xlabel('Real beta', fontsize=12)
    ax_beta.set_ylabel('Estimate of beta', fontsize=12)
    ax_beta.legend(loc="upper left", fontsize=8)



    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%stest_mean_estimate.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()



def predicted_observed_prevalence_truncated_gamma(min_minor_allele_cov=10):

    iter = 100

    x_mean_all = numpy.logspace(-2, -0.3, num=100)
    #x_mean_all  = [0.3]

    n_hosts = 100
    coverage = 80
    beta = 1

    x_mean_estimate_all = []
    observed_prevalence_all = []
    predicted_prevalence_all = []

    for x_mean_i in x_mean_all:

        shape = beta
        rate = beta/x_mean_i
        scale = 1/rate

        x_mean_estimate_all_i = []
        observed_prevalence_all_i = []
        predicted_prevalence_all_i = []
        for i in range(iter):

            # 50 fold coverage
            #D = numpy.repeat(50, n_hosts)
            D = numpy.random.poisson(lam=coverage, size=n_hosts)

            x = stats.gamma.rvs(shape, scale=scale, size=n_hosts)
            #x[x>1] = 1
            x_to_keep = x[x<=1]
            D_to_keep = D[x<=1]

            A = numpy.random.binomial(D_to_keep, x_to_keep)
            A[A < 10] = 0

            observed_prevalence = sum(A>0)/len(A)

            if observed_prevalence == 0:
                continue

            freqs = A/D_to_keep
            f_mean = numpy.mean(freqs)
            f_beta = (numpy.mean(freqs)**2)/numpy.var(freqs)

            prob_zero = (1+ ((f_mean/f_beta)*D_to_keep))**(-1*f_beta )

            if min_minor_allele_cov > 1:
                # range() excludes upper bound
                cumulative_sum = 0
                for n_i in range(0, min_minor_allele_cov):
                    cumulative_sum += (scipy.special.gamma(f_beta+n_i) / numpy.math.factorial(n_i)) * ((f_mean*D_to_keep/ (f_beta + f_mean*D_to_keep)) ** n_i) / scipy.special.gamma(f_beta)

                prob_zero *= cumulative_sum

            predicted_prevalence_slm = 1 - numpy.mean(prob_zero)

            x_mean_estimate_all_i.append(f_mean)
            observed_prevalence_all_i.append(observed_prevalence)
            predicted_prevalence_all_i.append(predicted_prevalence_slm)

        x_mean_estimate_all.append(numpy.mean(x_mean_estimate_all_i))
        observed_prevalence_all.append(numpy.mean(observed_prevalence_all_i))
        predicted_prevalence_all.append(numpy.mean(predicted_prevalence_all_i))


    x_mean_estimate_all = numpy.asarray(x_mean_estimate_all)
    observed_prevalence_all = numpy.asarray(observed_prevalence_all)
    predicted_prevalence_all = numpy.asarray(predicted_prevalence_all)

    relative_error = numpy.absolute(observed_prevalence_all - predicted_prevalence_all)/observed_prevalence_all


    fig = plt.figure(figsize = (8, 4)) #
    fig.subplots_adjust(bottom= 0.15)

    ax_prediction = plt.subplot2grid((1, 2), (0, 0), colspan=1)

    ax_error = plt.subplot2grid((1, 2), (0, 1), colspan=1)


    ax_prediction.scatter(observed_prevalence_all, predicted_prevalence_all)
    ax_prediction.set_xscale('log', basex=10)
    ax_prediction.set_yscale('log', basey=10)

    all_ = numpy.concatenate([observed_prevalence_all,predicted_prevalence_all])

    max_ = max(all_)*1.1
    min_ = min(all_)*0.8

    #print(min(f_max_no_zeros))

    ax_prediction.plot([min_, max_],[min_, max_], ls='--', lw=2, c='k', zorder=2)
    ax_prediction.set_xlim([min_, max_])
    ax_prediction.set_ylim([min_, max_])

    ax_prediction.set_xlabel('Observed prevalence', fontsize=12)
    ax_prediction.set_ylabel('Predicted prevalence', fontsize=12)

    ax_prediction.xaxis.set_tick_params(labelsize=7)
    ax_prediction.yaxis.set_tick_params(labelsize=7)



    ax_error.scatter(x_mean_estimate_all, relative_error)
    ax_error.set_xscale('log', basex=10)
    ax_error.set_yscale('log', basey=10)


    ax_error.set_xlim([min(x_mean_estimate_all), max(x_mean_estimate_all)])
    #ax_error.set_ylim([min_, max_])

    ax_error.set_xlabel('Observed mean frequency', fontsize=12)
    ax_error.set_ylabel('Mean relative error', fontsize=12)

    ax_error.xaxis.set_tick_params(labelsize=7)
    ax_error.yaxis.set_tick_params(labelsize=7)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.22, hspace=0.25)
    fig.savefig("%stest_gamma_truncation.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()




def test_numerical_estimate_slm(error_ = 0.01):

    n_hosts = 100
    coverage = 100

    D = numpy.random.poisson(lam=coverage, size=n_hosts)

    x_mean_observed = 0.01
    beta_observed = 0.9

    #rate = beta/x_mean_i
    #scale = 1/rate

    x_mean_0 = 0.04
    beta_0 = 1

    x_0 = stats.gamma.rvs(beta_0, scale=x_mean_0/beta_0, size=n_hosts)

    # truncate frequencies greater than 1
    x_to_keep = x_0[x_0<=1]
    D_to_keep = D[x_0<=1]

    A = numpy.random.binomial(D_to_keep, x_to_keep)
    A[A < 10] = 0
    freqs_truncated_0 = A/D_to_keep
    x_mean_truncated_0 = numpy.mean(freqs_truncated_0)
    beta_truncated_0 = (x_mean_truncated_0**2) / numpy.var(freqs_truncated_0)


    error_mean_x = numpy.absolute(x_mean_truncated_0 - x_mean_observed)/x_mean_observed
    error_beta = numpy.absolute(beta_truncated_0 - beta_observed)/beta_observed


    print(x_mean_truncated_0, x_mean_observed)
    print(beta_truncated_0, beta_observed)



    #print(error_mean_x, error_beta)



def test_joint_parameter_simulation(x_mean_observed, beta_observed, D, iter = 1000):

    #x_mean_all = numpy.logspace(-2, -0.4 , num=100)
    #beta_all = numpy.logspace(-2, 1 , num=100)

    n_hosts = len(D)

    x_mean_observed_log10 = numpy.log10(x_mean_observed)
    beta_observed_log10 = numpy.log10(beta_observed)

    x_mean_all_log10 = numpy.random.uniform(low=-4, high=-0.3, size=iter)
    beta_all_log10 = numpy.random.uniform(low=-2.3, high=1.5, size=iter)

    x_mean_all = 10**x_mean_all_log10
    beta_all = 10**beta_all_log10

    x_mean_truncated_all = []
    beta_truncated_all = []

    x_mean_all_to_keep = []
    beta_all_all_to_keep = []

    for i in range(iter):

        x_mean_i = x_mean_all[i]
        beta_i = beta_all[i]

        x_i = stats.gamma.rvs(beta_i, scale=x_mean_i/beta_i, size=n_hosts)
        x_i_to_keep = x_i[x_i<=1]
        D_i_to_keep = D[x_i<=1]
        A_i = numpy.random.binomial(D_i_to_keep, x_i_to_keep)
        A_i[A_i<10] = 0

        freqs_truncated_i = A_i/D_i_to_keep

        x_mean_truncated_i = numpy.mean(freqs_truncated_i)
        if x_mean_truncated_i == 0:
            continue
        beta_truncated_i = (x_mean_truncated_i**2) / numpy.var(freqs_truncated_i)

        x_mean_truncated_all.append(x_mean_truncated_i)
        beta_truncated_all.append(beta_truncated_i)

        x_mean_all_to_keep.append(x_mean_i)
        beta_all_all_to_keep.append(beta_i)


    x_mean_truncated_all = numpy.absolute(x_mean_truncated_all)
    beta_truncated_all = numpy.absolute(beta_truncated_all)

    x_mean_truncated_all_log10 = numpy.log10(x_mean_truncated_all)
    beta_truncated_all_log10 = numpy.log10(beta_truncated_all)

    #x_mean_truncated_all = numpy.sqrt(((x_mean_truncated_all - x_mean_observed)**2) + ((beta_truncated_all - beta_observed)**2))
    distance_log10 = numpy.sqrt(((x_mean_truncated_all_log10 - x_mean_observed_log10)**2) + ((beta_truncated_all_log10 - beta_observed_log10)**2))

    min_distance = numpy.amin(distance_log10)

    idx_ = numpy.where(distance_log10 == min_distance)[0][0]

    x_mean_truncated_best = x_mean_truncated_all[idx_]
    beta_truncated_best = beta_truncated_all[idx_]

    real_x_mean = x_mean_all_to_keep[idx_]
    real_beta = beta_all_all_to_keep[idx_]

    return x_mean_truncated_best, beta_truncated_best, min_distance







#n_hosts = 100
#coverage = 100
#D = numpy.random.poisson(lam=coverage, size=n_hosts)

#x_mean_observed = 0.01
#beta_observed = 0.05

#test_joint_parameter_simulation(x_mean_observed, beta_observed, D)
#print("hey")
#test_joint_parameter_simulation(x_mean_observed, beta_observed, D)
#test_joint_parameter_simulation(x_mean_observed, beta_observed, D)
#test_joint_parameter_simulation(x_mean_observed, beta_observed, D)

#test_gamma_parameter_sampling()

#test_numerical_estimate_slm()

# see if there's a single optimum
#def test_heatmap():



#test_numerical_estimate_slm()


#predicted_observed_prevalence_truncated_gamma()

#test_gamma_parameter_sampling()
