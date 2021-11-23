
import config
import numpy
import matplotlib.pyplot as plt



shape_ = numpy.random.uniform(low=1, high=10, size=100)
scale_ = numpy.random.uniform(low=1, high=10, size=100)
all_samples = []
for i in range(len(shape_)):

    shape_i = shape_[i]
    scale_i = scale_[i]

    sample = numpy.random.gamma(shape_i, scale=scale_i, size=100000)
    all_samples.append(sample)



all_samples_flat = numpy.concatenate(all_samples)
all_samples_flat = all_samples_flat[all_samples_flat>0]
all_samples_flat_log10 = numpy.log10(all_samples_flat)


fig, ax = plt.subplots(figsize=(4,4))


hist, bin_edges = numpy.histogram(all_samples_flat_log10, density=True, bins=20)
bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]

ax.scatter(bins_mean, hist, alpha=0.9, s=10)


fig.tight_layout()
fig.subplots_adjust(wspace=0.22, hspace=0.22)
fig.savefig("%sgamma_mix.png" % (config.analysis_directory), format='png', bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
plt.close()
