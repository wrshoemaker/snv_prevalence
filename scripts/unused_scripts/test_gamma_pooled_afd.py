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


import matplotlib.pyplot as plt
from scipy.stats import gamma



strain_directory = os.path.expanduser("~/GitHub/negative_selection_microbiome/data/strain_data_3/")

figure_directory = os.path.expanduser("~/GitHub/negative_selection_microbiome/analysis/")


afd_dict = {}

for filename in os.listdir(strain_directory):

    if filename.endswith(".pkl"):

        filepath = '%s/%s' % (strain_directory, filename)

        with open(filepath, 'rb') as handle:
            strain_dict = pickle.load(handle)

        #print(strain_dict)
        # exclue
        for species, freqs in strain_dict.items():

            # exclude samples with no strains
            if len(freqs) == 1:
                continue

            if species not in afd_dict:
                afd_dict[species] = []

            afd_dict[species].extend(freqs)




afd_all = [afd_dict[x] for x in afd_dict.keys() if len(afd_dict[x]) < 40]
afd_all = [item for sublist in afd_all for item in sublist]

mean_afd_all = numpy.mean(afd_all)
std_afd_all = numpy.std(afd_all)


fig, ax = plt.subplots(figsize=(4,4))

afd_all = []

count = 0
for species, afd in afd_dict.items():

    if len(afd) < 40:
        continue


    count += 1

    afd = numpy.asarray(afd)

    #afd = afd[afd>0.2]

    afd_log10 = numpy.log10(afd)

    rescaled_afd_log10 = (afd_log10 - mean_afd_all) / std_afd_all

    afd_all.extend(rescaled_afd_log10)


    bin_height, bin_boundary = numpy.histogram(rescaled_afd_log10, density=True, bins=10)

    #ax.hist(rescaled_afd_log10, alpha=0.8, bins= 20, density=True)
    # weights=numpy.zeros_like(rescaled_AFDs) + 1. / len(rescaled_AFDs)
    #ax.plot(x_range, gamma.pdf(x_range, ag, bg,cg), 'k', label='Gamma fit', lw=2)

    ax.scatter(bin_boundary[:-1], bin_height, alpha = 0.7)


#out = np.concatenate(input_list).ravel()

print(len(afd_all))


afd_all = numpy.asarray(afd_all)

ag,bg,cg = gamma.fit(10**afd_all)

x_range = numpy.linspace(min(10**afd_all) , max(10**afd_all) , 10000)
x_range_fit = gamma.pdf(x_range, ag, bg,cg)

#ax.plot(numpy.log10(x_range), x_range_fit, 'k', label='Gamma fit', lw=2)


ax.set_yscale('log', base=10)


ax.set_xlabel('Rescaled log\nrelative abundance', fontsize=12)
ax.set_ylabel('Probability density', fontsize=12)

ax.set_title("HMP strain data", fontsize=12)

fig.subplots_adjust(wspace=0.3, hspace=0.3)
fig.savefig(figure_directory + "afd.pdf", format='pdf', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
