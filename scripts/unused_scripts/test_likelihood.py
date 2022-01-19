from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import collections
import os.path
import scipy.stats as stats

import diversity_utils
import core_gene_utils
import parse_midas_data
import parse_HMP_data
import sample_utils
import calculate_substitution_rates
import clade_utils

import matplotlib.pyplot as plt
from scipy.stats import gamma
import scipy.special
from scipy.integrate import quad



def likelihood_poly(error_estimate, freq_estimate):

    
