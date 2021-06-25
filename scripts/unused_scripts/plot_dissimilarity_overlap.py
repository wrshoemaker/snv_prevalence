from __future__ import division
import os, sys
import bz2
import random
import itertools
import config
import parse_midas_data
import numpy
import pickle

from math import ceil
from scipy import linalg

import gzip

import matplotlib.pyplot as plt

from statsmodels import nonparametric
import statsmodels

import statsmodels.api as sm

data_directory = config.data_directory

species_name = 'Bacteroides_vulgatus_57955'


filename = "%sdissimilarity_overlap/%s.txt.gz" % (data_directory, species_name)





def lowess(x, y, f=2. / 3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(ceil(f * n))
    h = [numpy.sort(numpy.abs(x - x[i]))[r] for i in range(n)]
    w = numpy.clip(numpy.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yest = numpy.zeros(n)
    delta = numpy.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = numpy.array([numpy.sum(weights * y), numpy.sum(weights * y * x)])
            A = numpy.array([[numpy.sum(weights), numpy.sum(weights * x)],
                          [numpy.sum(weights * x), numpy.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]

        residuals = y - yest
        s = numpy.median(numpy.abs(residuals))
        delta = numpy.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2

    return yest






f = gzip.open(filename,'rb')


header = f.readline()

overlap_4D = []
dissimilarity_4D = []

overlap_1D = []
dissimilarity_1D = []
for line in f:

    line_split = line.strip().split(', ')
    if line_split[0] == '4D':

        overlap_4D.append(float(line_split[3]))
        dissimilarity_4D.append(float(line_split[4]))

    elif line_split[0] == '1D':
        overlap_1D.append(float(line_split[3]))
        dissimilarity_1D.append(float(line_split[4]))


    else:
        continue

f.close()


def calculate_bin_mean(overlap_, dissimilarity_):

    mean = []
    bin_midpoint = []

    hist, bin_edges = numpy.histogram(overlap_, bins=30)

    for bin_edge_idx in range(len(bin_edges)-1):

        lower_bin = bin_edges[bin_edge_idx]
        upper_bin = bin_edges[bin_edge_idx+1]

        dissimilarity_bin = dissimilarity_[ (overlap_ > lower_bin) & (overlap_ <= upper_bin) ]

        mean.append(numpy.mean(dissimilarity_bin))
        bin_midpoint.append((upper_bin+lower_bin)/2)

    return bin_midpoint, mean


overlap_4D = numpy.asarray(overlap_4D)
dissimilarity_4D = numpy.asarray(dissimilarity_4D)

overlap_1D = numpy.asarray(overlap_1D)
dissimilarity_1D = numpy.asarray(dissimilarity_1D)

bin_midpoint_4D, mean_4D = calculate_bin_mean(overlap_4D, dissimilarity_4D)

bin_midpoint_1D, mean_1D = calculate_bin_mean(overlap_1D, dissimilarity_1D)

f = 0.1
#yest_4D = lowess(overlap_4D, dissimilarity_4D, f=f, iter=3)
#yest_1D = lowess(overlap_1D, dissimilarity_1D, f=f, iter=3)


# bin y axis


fig, ax = plt.subplots(figsize = (4, 4))




ax.set_ylim([-0.05, 0.18])

#ax.scatter(overlap_4D, dissimilarity_4D, alpha = 0.05, color = 'b')
#ax.scatter(overlap_1D, dissimilarity_1D, alpha = 0.05, color = 'b')

ax.plot(bin_midpoint_4D, mean_4D, color = 'b')
ax.plot(bin_midpoint_1D, mean_1D, color = 'r')


fig.subplots_adjust(hspace=0.4, wspace=0.35) #hspace=0.3, wspace=0.5
fig_name =  "%sdissimilarity_overlap.png" % parse_midas_data.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
