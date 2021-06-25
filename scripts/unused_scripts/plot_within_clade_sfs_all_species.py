# SFS plots, comparing syn and non syn
# (a) regular SFS (percentage of population)
# (b) Singletons, doubletons, and pi-weighted

# Summary across species?
# dN/dS vs dS plot?
# core vs variable genes?

# Within-snp gene changes

import matplotlib
matplotlib.use('Agg')
import config
import parse_midas_data
###
#
# For today while the new data processes
#
import os
#parse_midas_data.data_directory = os.path.expanduser("~/ben_nandita_hmp_data_062517/")
#########################################
import parse_HMP_data


import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint


def get_sfs_diffusion_selection(f, gamma, theta=0.01):

    return theta * (1 / (f*(1-f))) * (1-exp(-2*gamma* (1-f) )) / (1-exp(-2*gamma))

def get_sfs_diffusion_neutral(f, theta=0.01):

    return theta/f


def get_sfs_diffusion_selection_folded(f, gamma, theta=0.01):

    return (theta * (1 / (f*(1-f))) * (1-exp(-2*gamma*f)) / (1-exp(-2*gamma))) + (theta * (1 / (f*(1-f))) * (1-exp(-2*gamma* (1-f) )) / (1-exp(-2*gamma)))

def get_sfs_diffusion_neutral_folded(f, theta=0.01):

    return theta/(f * (1-f))



def get_frequency(f, proportion_d, N_s_d):
    return (1-proportion_d) + (proportion_d * (exp(-1*N_s_d*f))  + exp(-1*N_s_d*(1-f)) )




def get_frequency_b_d(f, N,proportion_b, proportion_d, N_s_d):
    return (1-proportion_d-proportion_b) + (proportion_d * (exp(-1*N_s_d*f) + exp(-1*N_s_d*(1-f))))  + proportion_b


def get_frequency_intercept(f, proportion_d, N_s_d, intercept):
    return intercept + (1-proportion_d) + (proportion_d * (exp(-1*N_s_d*f))  + exp(-1*N_s_d*(1-f)) )





file = open("%smaf_all_species.txt" % (config.data_directory), 'r')

header_line = file.readline() # header
header_items = header_line.split(", ")

fig, ax = plt.subplots(figsize=(4,4))

ax.axhline(y=1, ls=':', color='grey')

for line in file:

    line = line.strip().split(', ')

    mafs =  numpy.asarray([float(x) for x in line[2].split(':')])

    sfs_synonymous = numpy.asarray([float(x) for x in line[4].split(':')])

    sfs_nonsynonymous = numpy.asarray([float(x) for x in line[3].split(':')])

    ax.plot(mafs, sfs_nonsynonymous/sfs_synonymous, linestyle='--', c='dodgerblue', alpha=0.6)


file.close()


ax.set_xlabel('Minor allele freq in largest clade, $f$')
ax.set_ylabel('Ratio of nonsynonymous and synonymous\nfraction of sites, $P(f_{N})/P(f_{S})$')
#ax.set_xlim([0,0.5])

ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)


#gamma = -10

#gammas = [-10, -1, -0.1]
#line_styles = ['-','--',':']
#for gamma_idx, gamma in enumerate(gammas):

#    sfs_ratio = [get_sfs_diffusion_selection_folded(maf_i, gamma=gamma) / get_sfs_diffusion_neutral_folded(maf_i) for maf_i in maf_range]
#    ax.plot(maf_range, sfs_ratio, linestyle=line_styles[gamma_idx] , c='k', label=r'$Ns=$' + str(gamma))

#ax.legend(loc="upper right", fontsize=8)

#ax.set_title( "All species", fontsize=12)#,y=0.95)
maf_range = numpy.linspace(0.01, 0.5, 1000)
sfs_ratio = [get_frequency(maf_i, 0.5, 10) for maf_i in maf_range]


#print(get_frequency_b_d(0, 100000, 0.1, 0.5, 10))


ax.plot(maf_range, sfs_ratio, linestyle='--' , c='k')#, label=r'$s_{b}/sp_{d}=%f$' % round(Ns_b/Ns_d, 2))


fig.savefig("%s%s.png" % (config.analysis_directory, 'within_clade_sfs_all_species'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



##### plot 2 s SFS ratio


#def get_frequency(f, proportion_b, proportion_d, N_s_b, N_s_d):
#    return (1-proportion_b-proportion_d) + (proportion_b * (exp(-1*N_s_b*f) + exp(-1*N_s_b*(1-f)))) + (proportion_d * (exp(-1*N_s_d*f) + exp(-1*N_s_d*(1-f)) ))
#    #denominator = (1-proportion_b-proportion_d)
#    #return 1 + (numerator / denominator)



#fig, ax = plt.subplots(figsize=(4,4))

#ax.axhline(y=1, ls=':', color='grey')

#maf_range = numpy.linspace(0.01, 0.5, 1000)

#population_scaled_selection = [0.01, 0.1, 1]
#line_styles = [':','--','-']

#for Ns_b in population_scaled_selection:
#    for Ns_d in population_scaled_selection:
#        Ns_d = -1* Ns_d
#Ns_b =  0
#Ns_d =  10
#sfs_ratio = [get_frequency(maf_i, 0.3, Ns_d) for maf_i in maf_range]

#if Ns_b > Ns_d:
#    color='b'
#elif Ns_b < Ns_d:
#    color='r'
#else:
#    color='darkgrey'


#ax.plot(maf_range, sfs_ratio, linestyle='--' , c=color)#, label=r'$s_{b}/sp_{d}=%f$' % round(Ns_b/Ns_d, 2))

#for gamma_idx, gamma in enumerate(gammas):

#ax.set_xlabel('Minor allele freq in largest clade, $f$')
#ax.set_ylabel('Ratio of nonsynonymous and synonymous\nfraction of sites, $P(f_{N})/P(f_{S})$')

#fig.savefig("%s%s.png" % (config.analysis_directory, 'sfs_ratio_model'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.close()




#sfs_axis.plot(mafs, synonymous_sfs*mafs*(1-mafs)/(synonymous_sfs*mafs*(1-mafs)).sum(), 'b.-',label='4D')
