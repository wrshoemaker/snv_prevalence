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

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
import scipy.stats as stats

x = numpy.random.normal(0.1, 1, size=100)

x = x[(x>0) & (x<1)]

b0 = 0.001
b1 = 2
y = b0*(x**b1)


print(x/y)
