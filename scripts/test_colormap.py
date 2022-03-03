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
from matplotlib import colors

import prevalence_utils

from collections import Counter

numpy.random.seed(123456789)


species_to_run = prevalence_utils.species_to_run
genera = [s.split('_')[0] for s in species_to_run]
count_genera = Counter(genera)



#color_map = {'Bacteroides': 'Blues', 'Alistipes': 'Greens', 'Ruminococcus': 'Oranges', ''}

orange_rgb = colors.hex2color(colors.cnames['orange'])

print(orange_rgb)

species_to_run = {'Alistipes_finegoldii_56071': colors.hex2color(colors.cnames['lightgreen']),
                    'Alistipes_onderdonkii_55464': colors.hex2color(colors.cnames['mediumseagreen']),
                    'Alistipes_putredinis_61533': colors.hex2color(colors.cnames['seagreen']),
                    'Alistipes_shahii_62199':  colors.hex2color(colors.cnames['darkgreen']),
                    'Bacteroidales_bacterium_58650': colors.hex2color(colors.cnames['khaki']),
                    'Bacteroides_caccae_53434': colors.hex2color(colors.cnames['skyblue']),
                    'Bacteroides_cellulosilyticus_58046': colors.hex2color(colors.cnames['aquamarine']),
                    'Bacteroides_fragilis_54507': colors.hex2color(colors.cnames['steelblue']),
                    'Bacteroides_ovatus_58035': colors.hex2color(colors.cnames['cyan']),
                    'Bacteroides_stercoris_56735': colors.hex2color(colors.cnames['lightskyblue']),
                    'Bacteroides_thetaiotaomicron_56941': colors.hex2color(colors.cnames['darkturquoise']),
                    'Bacteroides_uniformis_57318': colors.hex2color(colors.cnames['teal']),
                    'Bacteroides_vulgatus_57955': colors.hex2color(colors.cnames['dodgerblue']),
                    'Bacteroides_xylanisolvens_57185': colors.hex2color(colors.cnames['cadetblue']),
                    'Barnesiella_intestinihominis_62208': colors.hex2color(colors.cnames['darkorange']),
                    'Dialister_invisus_61905': colors.hex2color(colors.cnames['sandybrown']),
                    'Eubacterium_rectale_56927': colors.hex2color(colors.cnames['olive']),
                    'Oscillibacter_sp_60799': colors.hex2color(colors.cnames['saddlebrown']),
                    'Parabacteroides_distasonis_56985': colors.hex2color(colors.cnames['orchid']),
                    'Parabacteroides_merdae_56972': colors.hex2color(colors.cnames['darkmagenta']),
                    'Ruminococcus_bicirculans_59300': colors.hex2color(colors.cnames['orangered']),
                    'Ruminococcus_bromii_62047': colors.hex2color(colors.cnames['darkred'])}

print(species_to_run)


#print(dict(count_genera))
