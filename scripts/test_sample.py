from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path
import random
from collections import Counter

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import plot_utils
import prevalence_utils

from scipy.stats import gamma, gaussian_kde
import scipy.stats as stats

import calculate_predicted_prevalence_mapgd



species_name = 'Prevotella_copri_61740'


directory = config.data_directory+"strain_data/"

for filename in os.listdir(directory):

    if filename.endswith(".pkl"):
        path = os.path.join(directory, filename)

        with open(path, 'rb') as handle:
            b = pickle.load(handle)

        if species_name in b:

            print(filename, len(b[species_name]))
