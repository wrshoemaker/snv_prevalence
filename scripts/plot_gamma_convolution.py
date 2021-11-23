from __future__ import division
import config
import numpy
import matplotlib.pyplot as plt

rel_abundances = [0.4, 0.2, 0.1, 0.05, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02, 0.1 ]

N = 1000000
mu = 10**-6
t = 10000

#x_1 = 0.2
#x_2 = 1-x_1

all_freqs = []
for r in rel_abundances:
    f_max = t / (N*r)
    theta = 2*mu*N*r
    freqs = numpy.random.gamma(theta, scale=f_max, size=1000)
    all_freqs.extend(freqs.tolist())

    print(f_max, theta)


all_freqs = numpy.asarray(all_freqs)
#freqs_log10 = numpy.log10(all_freqs)
#freqs_log10_rescaled = (freqs_log10 - numpy.mean(freqs_log10)) / numpy.std(freqs_log10)



#freqs_2 = numpy.random.gamma(theta_2, scale=f_max_2, size=1000)


#freqs = numpy.concatenate([freqs_1, freqs_2])



fig, ax = plt.subplots(figsize=(4,4))

hist, bin_edges = numpy.histogram(all_freqs, density=True, bins=20)
bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]

ax.scatter(bins_mean, hist, alpha=0.9, s=10)

fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.22)
fig.savefig("%sgamma_convolution.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
