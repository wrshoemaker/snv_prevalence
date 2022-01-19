from __future__ import division
import sys
import numpy
import pickle
import bz2
import gzip
import config
import math
import os.path

import diversity_utils
import figure_utils
import parse_midas_data
import prevalence_utils
import scipy.stats as stats
import scipy.special as special

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import calculate_predicted_prevalence_mapgd
import calculate_predicted_prevalence



def prob_zero_binomial(f,D):
    return (1-f)**D


def prob_zero_poisson(f,D):
    return numpy.exp(-1*D*f)


f = 0.01
D=50
#print(1-prob_zero_binomial(f,D))
#print(1-prob_zero_poisson(f,D))

# any error due to the Poisson approximation would cause an underestimation of prevalence.

beta = 0.8
x_mean = 0.6
N = 100


def prob_zero(depth):
    return (1 + (x_mean/beta)*depth)**(-1*beta)

def prob_zero_truncated():

    return prob_zero() / special.gammainc(beta , beta/x_mean)



#print(1-prob_zero())
#print(1-prob_zero_truncated())


print(1-prob_zero(30))

print(1-prob_zero(23))



# test jacopos approach


freqs = numpy.asarray([0.1, 0.05, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0])
depths = numpy.asarray([20, 50, 35, 23, 54, 54, 54, 54, 40, 52, 40, 36])

freqs_no_zeros = freqs[freqs>0]

tf = numpy.mean(freqs_no_zeros/ depths[freqs>0])
tvpf = numpy.mean( (freqs_no_zeros**2 - freqs_no_zeros) / (depths[freqs>0]**2) )

prevalence = sum(freqs>0)/len(freqs)

f = prevalence*tf
vf= prevalence*tvpf

# there's this command in Jacopo's code %>% mutate(vf = vf - f^2 )%>%
# It's applied after f and vf are calculated
# This should be equivalent to the mean and variance including zero
vf = vf - (f**2)

beta = (f**2)/vf
theta = f/beta


#print(theta, beta)
#predicted_prevalence_slm = 1 - numpy.mean(((1+theta*depths)**(-1*beta )))



#prob_absence = (1+theta*depths)**(-1*beta )




#print(prob_absence)

#print(special.gammainc(0.001 , 0))






x_range = numpy.linspace(-4, 3, 10000)

k = 2.0
k_digamma = special.digamma(k)
k_trigamma = special.polygamma(1,k)
gammalog = k*k_trigamma*x_range - numpy.exp(numpy.sqrt(k_trigamma)*x_range + k_digamma) - numpy.log(special.gamma(k)) + k*k_digamma + numpy.log10(numpy.exp(1))


fig, ax = plt.subplots(figsize=(4,4))
ax.plot(x_range, 10**gammalog, 'k', label='Gamma', lw=2)
ax.set_ylim([0.007, 1])
ax.set_yscale('log', basey=10)


fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.25)
fig.savefig("%sgamma_dist.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
