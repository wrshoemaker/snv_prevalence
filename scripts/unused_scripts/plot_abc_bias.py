from __future__ import division
import os
import sys
import json
import pickle
import numpy
import config
import itertools

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import scipy.stats as stats

from scipy.stats import gamma, gaussian_kde

from sklearn import svm
from sklearn.model_selection import train_test_split


#import simulate_gamma_estimator

min_f_mean = 2.304437869822499e-05
max_f_mean = 0.39
min_beta = 0.00278551532033426
max_beta = 5.874376635200855


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
    idx_error = (numpy.isfinite(error_ratio)) & (beta_truncated_all>min_beta) & (beta_truncated_all<max_beta)  & (f_mean_truncated_all>min_f_mean) & (f_mean_truncated_all<max_f_mean)

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




f_mean_truncated_all, beta_truncated_all, error_ratio, f_mean_true_all, beta_true_all, f_mean_abc_all, beta_abc_all, predicted_prevalence_truncated_all, observed_prevalence_truncated_all, predicted_prevalence_abc_all, observed_prevalence_abc_all = get_simulated_data()

f_mean_bias_naive = (f_mean_truncated_all - f_mean_true_all )
beta_bias_naive = (beta_truncated_all - beta_true_all)

f_mean_bias_abc = (f_mean_abc_all - f_mean_true_all)
beta_bias_abc = (beta_abc_all - beta_true_all)

all_f_mean_biases = list(itertools.chain.from_iterable( [f_mean_bias_naive,  f_mean_bias_abc] ))
all_beta_biases = list(itertools.chain.from_iterable( [beta_bias_naive,  beta_bias_abc] ))


error_f_mean_naive = numpy.absolute(f_mean_truncated_all - f_mean_true_all ) / f_mean_true_all
error_beta_naive = numpy.absolute(beta_truncated_all - beta_true_all ) / beta_true_all


#print(f_mean_truncated_all[:10], f_mean_true_all[:10])


# heatmaps

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar






def run_svm(x, y, z):

    log_z = numpy.log10(z)

    # run SVM
    dep_variable = numpy.copy(log_z)
    #y[y<=0] = 0
    #y[y>0] = 1

    dep_variable[dep_variable>=0] = False
    dep_variable[dep_variable<0] = True

    # 0 = dont run ABC
    # 1 = run ABC

    x_log10 = numpy.log10(x)
    y_log10 = numpy.log10(y)

    indep_variable = numpy.stack([x_log10, y_log10], axis=1)

    # select random observations
    rndm_idx = numpy.random.choice(len(dep_variable), size=10000, replace=False)

    indep_rndm = indep_variable[rndm_idx,:]
    dep_rndm = dep_variable[rndm_idx]

    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)

    clf = svm.SVC(kernel='rbf', gamma='auto', C=1)
    #clf = svm.SVC(kernel='poly', gamma='auto', C=0.1)
    clf.fit(indep_rndm, dep_rndm)

    return clf








def plot_estimator_bias_abc():

    gs = gridspec.GridSpec(nrows=2, ncols=2)
    fig = plt.figure(figsize = (9, 8))
    ax_f_mean_naive = fig.add_subplot(gs[0, 0])
    ax_beta_naive = fig.add_subplot(gs[0, 1])
    ax_f_mean_abc = fig.add_subplot(gs[1, 0])
    ax_beta_abc = fig.add_subplot(gs[1, 1])

    axes = [ax_f_mean_naive, ax_beta_naive, ax_f_mean_abc, ax_beta_abc]
    bias_all = [f_mean_bias_naive, beta_bias_naive, error_f_mean_naive, error_beta_naive]


    cbar_labels = ['Bias of naive estimate of ' + r'$\bar{f}$', 'Bias of naive estimate of ' + r'$\beta$', 'Relative error of naive estimate of ' + r'$\bar{f}$' + ', ' + r'$\epsilon$', 'Relative error of naive estimate of ' + r'$\beta$' +  ', ' + r'$\epsilon$']

    for bias_idx, bias in enumerate(bias_all):

        ax = axes[bias_idx]

        #if (bias_idx == 0) or (bias_idx == 2):
        #    min_cbar = min(all_f_mean_biases)
        #    max_cbar = max(all_f_mean_biases)

        min_ =  max(numpy.absolute(bias)) *-1
        max_ =  max(numpy.absolute(bias))

        #else:
        #    min_cbar = min(all_beta_biases)
        #    max_cbar = max(all_beta_biases)

        if (bias_idx == 0) or (bias_idx == 1):

            s = ax.scatter(f_mean_true_all, beta_true_all, edgecolors='none', s=10, alpha=0.6, c=bias, norm=colors.Normalize(vmin=min_, vmax=max_), cmap='RdBu_r')

        else:
            f_mean_true_all_to_plot = f_mean_true_all[bias>0]
            beta_true_all_to_plot = beta_true_all[bias>0]
            bias = bias[bias>0]
            s = ax.scatter(f_mean_true_all_to_plot, beta_true_all_to_plot, edgecolors='none', s=10, alpha=1, c=bias, norm=colors.LogNorm(vmin=min(bias), vmax=max(bias)), cmap='Blues')

            #s = ax.scatter(f_mean_truncated_all_to_plot, beta_truncated_all_to_plot, edgecolors='none', s=numpy.absolute(log_distance_ratio), alpha=0.6, c=distance_ratio, norm=colors.LogNorm(), cmap='RdBu_r')

        cbar = plt.colorbar(mappable=s, ax=ax)
        cbar.set_label(cbar_labels[bias_idx])

        ax.set_xlabel('True mean frequency, ' + r'$\bar{f}$', fontsize=11)
        ax.set_ylabel('True squared inverse CV, ' + r'$\beta$', fontsize=11)

        ax.set_xlim([min(f_mean_true_all), max(f_mean_true_all)])
        ax.set_ylim([min(beta_true_all), max(beta_true_all)])

        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%sgamma_estimator_bias.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
    plt.close()




def plot_estimator_comparison():

    idx = (f_mean_truncated_all > 0) & (f_mean_true_all>0) & (beta_truncated_all>0) & (beta_true_all>0) & (f_mean_abc_all>0) & (beta_abc_all>0)

    f_mean_truncated_all_to_keep = f_mean_truncated_all[idx]
    beta_truncated_all_to_keep = beta_truncated_all[idx]

    f_mean_true_all_to_keep = f_mean_true_all[idx]
    beta_true_all_to_keep = beta_true_all[idx]

    f_mean_abc_all_to_keep = f_mean_abc_all[idx]
    beta_abc_all_to_keep = beta_abc_all[idx]

    distance_naive = numpy.sqrt(  ((numpy.log10(f_mean_truncated_all_to_keep) - numpy.log10(f_mean_true_all_to_keep))**2) + (numpy.log10(beta_truncated_all_to_keep) - numpy.log10(beta_true_all_to_keep))**2   )
    distance_abc = numpy.sqrt(  ((numpy.log10(f_mean_abc_all_to_keep) - numpy.log10(f_mean_true_all_to_keep))**2) + (numpy.log10(beta_abc_all_to_keep) - numpy.log10(beta_true_all_to_keep))**2   )

    error_f_mean_naive = numpy.absolute(f_mean_truncated_all_to_keep - f_mean_true_all_to_keep)/f_mean_true_all_to_keep
    error_beta_naive = numpy.absolute(beta_truncated_all_to_keep - beta_true_all_to_keep)/beta_truncated_all_to_keep

    error_f_mean_abc = numpy.absolute(f_mean_abc_all_to_keep - f_mean_true_all_to_keep)/f_mean_true_all_to_keep
    error_beta_abc = numpy.absolute(beta_abc_all_to_keep - beta_true_all_to_keep)/beta_truncated_all_to_keep

    error_sum_naive = error_f_mean_naive + error_beta_naive
    error_sum_abc = error_f_mean_abc + error_beta_abc

    #idx = (~numpy.isnan(distance_naive)) & (~numpy.isnan(distance_abc)) & (distance_naive > 0) & (distance_abc > 0)

    #idx = (error_sum_naive > 0) & (error_sum_abc > 0)
    idx = (distance_naive > 0) & (distance_abc > 0)

    distance_naive = distance_naive[idx]
    distance_abc = distance_abc[idx]

    f_mean_truncated_all_to_plot = f_mean_truncated_all_to_keep[idx]
    beta_truncated_all_to_plot = beta_truncated_all_to_keep[idx]

    distance_ratio = distance_abc/distance_naive

    log_distance_ratio = numpy.log10(distance_ratio)


    fig, ax = plt.subplots(figsize=(5,4))

    s = ax.scatter(f_mean_truncated_all_to_plot, beta_truncated_all_to_plot, edgecolors='none', s=numpy.absolute(log_distance_ratio), alpha=0.6, c=distance_ratio, norm=colors.LogNorm(), cmap='RdBu_r')

    cbar = plt.colorbar(mappable=s, ax=ax)
    cbar.set_label('Ratio of distances, ' +  r'$d_{\mathrm{ABC}}/d_{\mathrm{Naive}}$')

    ax.set_xlabel('True mean frequency, ' + r'$\bar{f}$', fontsize=11)
    ax.set_ylabel('True squared inverse CV, ' + r'$\beta$', fontsize=11)

    #ax.set_xlim([min(f_mean_truncated_all_to_plot), max(f_mean_truncated_all_to_plot)])
    #ax.set_ylim([min(beta_truncated_all_to_plot), max(beta_truncated_all_to_plot)])

    f_mean_truncated_all_to_plot_log10 = numpy.log10(f_mean_truncated_all_to_plot)
    beta_truncated_all_to_plot_log10 = numpy.log10(beta_truncated_all_to_plot)

    clf = run_svm(f_mean_truncated_all_to_plot, beta_truncated_all_to_plot, distance_ratio)


    XX, YY = numpy.mgrid[min(f_mean_truncated_all_to_plot_log10):max(f_mean_truncated_all_to_plot_log10):200j, min(beta_truncated_all_to_plot_log10):max(beta_truncated_all_to_plot_log10):200j]
    Z = clf.decision_function(numpy.c_[XX.ravel(), YY.ravel()])
    #Put the result into a color plot
    Z = Z.reshape(XX.shape)

    ax.contour(10**XX, 10**YY, Z, colors=["k", "k", "k"], linestyles=["--", "-", "--"], levels=[-0.5, 0, 0.5])
    #ax.contour(10**XX, 10**YY, Z, colors=["k"], linestyles=["-"], levels=[0])

    # '--': support vector
    # '-': separating hyperplane

    legend_elements = [Line2D([0], [0], color='k', lw=2, ls='-', label='Separating hyperplane'),
                       Line2D([0], [0], color='k', lw=2, ls='--', label='Support vector')]

    ax.legend(loc="upper left", handles=legend_elements, fontsize=8)

    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%sestimator_distance_comparison.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
    plt.close()







def plot_abc_prevalence():

    gs = gridspec.GridSpec(nrows=1, ncols=3)
    fig = plt.figure(figsize = (12, 4))
    ax_naive = fig.add_subplot(gs[0, 0])
    ax_abc =fig.add_subplot(gs[0, 1])
    ax_vs =fig.add_subplot(gs[0, 2])

    f_mean_truncated_all, beta_truncated_all, error_ratio, f_mean_true_all, beta_true_all, f_mean_abc_all, beta_abc_all, predicted_prevalence_truncated_all, observed_prevalence_truncated_all, predicted_prevalence_abc_all, observed_prevalence_abc_all = get_simulated_data()



    idx_to_plot = numpy.random.choice(len(observed_prevalence_truncated_all), size=10000, replace=False)

    observed_prevalence_truncated_all = observed_prevalence_truncated_all[idx_to_plot]
    predicted_prevalence_truncated_all = predicted_prevalence_truncated_all[idx_to_plot]
    observed_prevalence_abc_all = observed_prevalence_abc_all[idx_to_plot]
    predicted_prevalence_abc_all = predicted_prevalence_abc_all[idx_to_plot]

    error_naive = numpy.absolute(observed_prevalence_truncated_all - predicted_prevalence_truncated_all)/observed_prevalence_truncated_all
    error_abc = numpy.absolute(observed_prevalence_abc_all - predicted_prevalence_abc_all)/observed_prevalence_abc_all

    idx_ = (error_abc>0) & (error_naive>0)

    observed_prevalence_truncated_all = observed_prevalence_truncated_all[idx_]
    predicted_prevalence_truncated_all = predicted_prevalence_truncated_all[idx_]
    observed_prevalence_abc_all = observed_prevalence_abc_all[idx_]
    predicted_prevalence_abc_all = predicted_prevalence_abc_all[idx_]

    error_naive = error_naive[idx_]
    error_abc = error_abc[idx_]

    #all_ = numpy.concatenate([predicted_prevalence_no_zeros, observed_prevalence_no_zeros])
    xy = numpy.vstack([observed_prevalence_truncated_all, predicted_prevalence_truncated_all])
    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = observed_prevalence_truncated_all[idx], predicted_prevalence_truncated_all[idx], z[idx]
    ax_naive.scatter(x, y, c=z, cmap="Blues", s=90, alpha=0.9, edgecolor='', zorder=1)

    #ax_naive.scatter(observed_prevalence_truncated_all, predicted_prevalence_truncated_all, alpha=0.2)

    min_naive = min([ min(observed_prevalence_truncated_all), min(predicted_prevalence_truncated_all) ])
    max_naive = max([ max(observed_prevalence_truncated_all), max(predicted_prevalence_truncated_all) ])

    ax_naive.plot([min_naive, max_naive],[min_naive, max_naive], 'k--')

    ax_naive.set_xscale('log', basex=10)
    ax_naive.set_yscale('log', basey=10)



    xy = numpy.vstack([observed_prevalence_abc_all, predicted_prevalence_abc_all])
    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = observed_prevalence_abc_all[idx], predicted_prevalence_abc_all[idx], z[idx]
    ax_abc.scatter(x, y, c=z, cmap="Blues", s=90, alpha=0.9, edgecolor='', zorder=1)
    #ax_abc.scatter(observed_prevalence_abc_all, predicted_prevalence_abc_all, alpha=0.2)

    min_abc = min([ min(observed_prevalence_abc_all), min(predicted_prevalence_abc_all) ])
    max_abc = max([ max(observed_prevalence_abc_all), max(predicted_prevalence_abc_all) ])

    ax_abc.plot([min_abc, max_abc],[min_abc, max_abc], 'k--')

    ax_abc.set_xscale('log', basex=10)
    ax_abc.set_yscale('log', basey=10)



    xy = numpy.vstack([error_naive, error_abc])
    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = error_naive[idx], error_abc[idx], z[idx]
    ax_vs.scatter(x, y, c=z, cmap="Blues", s=90, alpha=0.9, edgecolor='', zorder=1)

    #ax_vs.scatter(observed_prevalence_abc_all, error_abc, alpha=0.2)

    ax_vs.set_xlim(min(error_naive), max(error_naive))
    ax_vs.set_ylim(min(error_abc), max(error_abc))

    #ax_vs.set_xlim(min(observed_prevalence_abc_all), max(observed_prevalence_abc_all))
    #ax_vs.set_ylim(min(error_abc), max(error_abc))


    #ax_vs.plot([min(error_naive), max(error_naive)],[min(error_abc), max(error_abc)], 'k--')

    ax_vs.set_xscale('log', basex=10)
    ax_vs.set_yscale('log', basey=10)


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2)
    # dpi = 600
    fig.savefig("%sabc_prevalence.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi=600)
    plt.close()





#plot_abc_prevalence()
plot_estimator_bias_abc()

#plot_estimator_comparison()

#def get_data(estimator='mean'):






#all_cdf = np.asarray(all_cdf)
#y_axis = np.transpose(y_axis)
