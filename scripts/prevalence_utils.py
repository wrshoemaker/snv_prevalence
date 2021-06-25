from __future__ import division
import sample_utils
import config
import parse_midas_data
import os.path
import os
import pylab
import sys
import numpy
import gzip
import pickle
import bz2

from scipy import stats
from scipy.stats import t

import matplotlib.cm as cm


conf=0.95


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']

good_species = 'Eubacterium_rectale_56927'
bad_species = 'Bacteroides_vulgatus_57955'

good_bad_color_dict = {good_species: 'dodgerblue', bad_species: '#FF6347'}
good_bad_color_map_dict = {good_species: cm.Blues, bad_species: cm.Reds}



def permutational_two_sample_t_test_equal_n(x, y, n=10000):


    def calculate_t(x, y):

        s_2_x = numpy.var(x, ddof=1)
        s_2_y = numpy.var(y, ddof=1)
        # pooled standard deviation
        s_p = numpy.sqrt((s_2_x+s_2_y)/2)

        return  (numpy.absolute(numpy.mean(x)) - numpy.absolute(numpy.mean(y))) / (numpy.sqrt(2/len(x)) * s_p)

    if type(x) is not numpy.ndarray:
        x = numpy.asarray(x)

    if type(y) is not numpy.ndarray:
        y = numpy.asarray(y)

    x_and_y = numpy.concatenate([x,y])

    t = calculate_t(x, y)

    t_abs_null = []

    for n_i in range(n):

        x_and_y_null = numpy.random.permutation(x_and_y)

        t_null_i = calculate_t(x_and_y_null[:len(x)], x_and_y_null[len(x):])

        t_abs_null.append(numpy.absolute(t_null_i))


    t_abs_null = numpy.asarray(t_abs_null)


    p = sum(t_abs_null > numpy.absolute(t)) / n

    return t, p











# performs regression on the two arrays passes
def get_confidence_hull(x, y):

    if type(x) is not numpy.ndarray:
        x = numpy.asarray(x)

    if type(y) is not numpy.ndarray:
        y = numpy.asarray(y)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)


    if min([min(x), min(y) ] ) < 0:
        min_range = min([min(x), min(y) ] ) * 1.3

    else:
        min_range = min([min(x), min(y) ] ) * 0.5

    max_range = max([max(x), max(y) ] ) * 1.3

    x_range = numpy.linspace(min_range, max_range, num=1000)
    y_range_pred = numpy.asarray([ intercept + (x_i*slope) for x_i in  x_range])

    y_pred = numpy.asarray([intercept + (slope*x_i) for x_i in x])

    SSE = sum((y - y_pred) ** 2)
    N = len(x)
    sd_SSE = numpy.sqrt( (1/ (N-2)) * SSE)
    sxd = numpy.sum((x-numpy.mean(x))**2)

    sx = (x_range-numpy.mean(x))**2	# x axisr for band
    # Quantile of Student's t distribution for p=1-alpha/2
    alpha = 1-conf
    q = stats.t.ppf(1-alpha/2, N-2)
    # Confidence band
    dy = q*sd_SSE*numpy.sqrt( 1/N + sx/sxd )
    # Upper confidence band
    ucb = y_range_pred + dy
    # Lower confidence band
    lcb = y_range_pred - dy


    return x_range, y_range_pred, lcb, ucb
