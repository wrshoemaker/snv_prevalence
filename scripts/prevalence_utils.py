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
import calculate_predicted_prevalence_mapgd

from scipy import stats
from scipy.stats import t

import matplotlib.cm as cm


conf=0.95

n_points=1000
color_radius=2


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']

good_species = 'Eubacterium_rectale_56927'
bad_species = 'Bacteroides_vulgatus_57955'

good_bad_color_dict = {good_species: 'dodgerblue', bad_species: '#FF6347'}
good_bad_color_map_dict = {good_species: cm.Blues, bad_species: cm.Reds}

variant_color_dict = {'4D': 'dodgerblue', '1D': '#FF6347'}

variant_cmap_dict = {'4D': 'Blues', '1D': 'Reds'}






species_to_run = ['Alistipes_finegoldii_56071', 'Alistipes_onderdonkii_55464', 'Alistipes_putredinis_61533',
                    'Alistipes_shahii_62199', 'Bacteroidales_bacterium_58650', 'Bacteroides_caccae_53434',
                    'Bacteroides_cellulosilyticus_58046', 'Bacteroides_fragilis_54507', 'Bacteroides_ovatus_58035',
                    'Bacteroides_stercoris_56735', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_uniformis_57318',
                    'Bacteroides_vulgatus_57955', 'Bacteroides_xylanisolvens_57185', 'Barnesiella_intestinihominis_62208',
                    'Dialister_invisus_61905', 'Eubacterium_rectale_56927', 'Oscillibacter_sp_60799', 'Parabacteroides_distasonis_56985',
                    'Parabacteroides_merdae_56972', 'Ruminococcus_bicirculans_59300', 'Ruminococcus_bromii_62047']




def get_relative_richness_dict(variant_type='4D'):

    max_richness = 4
    richness_range = list(range(1, max_richness+1))

    good_species_list = calculate_predicted_prevalence_mapgd.good_species_list

    relative_richness_dict = {}

    number_samples_dict = {}

    intermediate_strain_filename_template = config.data_directory+"strain_data/%s.pkl"

    for species_name in good_species_list:

        pi_dict = calculate_predicted_prevalence_mapgd.load_pi_dict(species_name)

        samples = list(pi_dict[variant_type].keys())

        number_samples_dict[species_name] = len(samples)

        richness_all = []

        for sample in samples:

            intermediate_strain_filename = intermediate_strain_filename_template % sample

            if os.path.isfile(intermediate_strain_filename) == False:
                continue

            with open(intermediate_strain_filename, 'rb') as handle:
                b = pickle.load(handle)

            if species_name in b:

                abundances = b[species_name]
                richness_all.append(len(abundances))


        relative_richness_dict[species_name] = {}

        for richness_i_idx, richness_i in enumerate(richness_range):

            relative_richness_dict[species_name][richness_i] = richness_all.count(richness_i)/len(richness_all)


    relative_richness_1 = [relative_richness_dict[s][1] for s in good_species_list]
    good_species_list_sorted = [s[0] for s in sorted(zip(good_species_list, relative_richness_1), key = lambda t: t[1])]#[::-1]
    proprtion_richness = [[],[],[],[]]

    for species_name in good_species_list_sorted:

        for richness_i_idx, richness_i in enumerate(richness_range):

            proprtion_richness[richness_i_idx].append(relative_richness_dict[species_name][richness_i])


    return relative_richness_dict, number_samples_dict, proprtion_richness






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





# https://github.com/weecology/macroecotools/blob/master/macroecotools/macroecotools.py
# code to cluster points
def count_pts_within_radius(x, y, radius, logscale=0):
    """Count the number of points within a fixed radius in 2D space"""
    #TODO: see if we can improve performance using KDTree.query_ball_point
    #http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_point.html
    #instead of doing the subset based on the circle
    unique_points = set([(x[i], y[i]) for i in range(len(x))])
    count_data = []
    logx, logy, logr = numpy.log10(x), numpy.log10(y), numpy.log10(radius)
    for a, b in unique_points:
        if logscale == 1:
            loga, logb = numpy.log10(a), numpy.log10(b)
            num_neighbors = len(x[((logx - loga) ** 2 +
                                   (logy - logb) ** 2) <= logr ** 2])
        else:
            num_neighbors = len(x[((x - a) ** 2 + (y - b) ** 2) <= radius ** 2])
        count_data.append((a, b, num_neighbors))
    return count_data


def plot_color_by_pt_dens(x, y, radius, loglog=0):
    """Plot bivariate relationships with large n using color for point density

    Inputs:
    x & y -- variables to be plotted
    radius -- the linear distance within which to count points as neighbors
    loglog -- a flag to indicate the use of a loglog plot (loglog = 1)

    The color of each point in the plot is determined by the logarithm (base 10)
    of the number of points that occur with a given radius of the focal point,
    with hotter colors indicating more points. The number of neighboring points
    is determined in linear space regardless of whether a loglog plot is
    presented.
    """
    plot_data = count_pts_within_radius(x, y, radius, loglog)
    sorted_plot_data = numpy.array(sorted(plot_data, key=lambda point: point[2]))

    return sorted_plot_data
