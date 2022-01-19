
from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config


output_filename = config.data_directory+"strain_diversity.dat"

with open(output_filename, 'rb') as handle:
    strain_diversity_dict = pickle.load(handle)


print(strain_diversity_dicts)
