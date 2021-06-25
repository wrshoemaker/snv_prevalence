import numpy
import matplotlib.pyplot as plt
import config


slope = 2
intercept = 1

x = numpy.random.uniform(low=0.0, high=1.0, size=10000)
y = slope*x + intercept

# increments of x
x_range = numpy.linspace(0, 1, num = 50)

y_from_max_xrange = []
max_xrange = []

y_from_min_xrange = []
min_xrange = []

for x_range_i in x_range:

    y_from_max_xrange_i = y[x > x_range_i]

    if len(y_from_max_xrange_i) > 2:

        y_from_max_xrange.append(numpy.mean(y_from_max_xrange_i))
        max_xrange.append(x_range_i)


    y_from_min_xrange_i = y[x < x_range_i]

    if len(y_from_min_xrange_i) > 2:

        y_from_min_xrange.append(numpy.mean(y_from_min_xrange_i))
        min_xrange.append(x_range_i)




fig, ax = plt.subplots(figsize=(4,4))

ax.scatter(max_xrange, y_from_max_xrange)




ax.set_xlabel('Minimum value of independent variable, ' + r'$x_{min}$', fontsize=10)
ax.set_ylabel('Mean truncated dependent variable, ' +  r'$y(x > x_{min})$', fontsize=9)


fig.savefig(config.analysis_directory + "truncat_lower.png", format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()




fig, ax = plt.subplots(figsize=(4,4))

ax.scatter(min_xrange, y_from_min_xrange)

ax.set_xlabel('Maximum value of independent variable, ' + r'$x_{max}$', fontsize=10)
ax.set_ylabel('Mean truncated dependent variable, ' +  r'$y(x < x_{max})$', fontsize=9)


fig.savefig(config.analysis_directory + "truncat_upper.png", format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()





fig, ax = plt.subplots(figsize=(4,4))

f0 = 0.05

original_weight = numpy.exp(-1 * x_range / f0)
modified_weight = numpy.exp(-1 * (1-x_range) / f0)

#ax.scatter(min_xrange, y_from_min_xrange)

ax.text(0.50,0.75, r'$f_{0}=0.05$', fontsize=9, color='k', ha='center', va='center', transform=ax.transAxes )


ax.plot(x_range, original_weight, c='b', label=r'$ e^{- \frac{f}{f_{0}}}$')
ax.plot(x_range, modified_weight, c='r', label=r'$ e^{- \frac{1-f}{f_{0}}}$')

ax.set_xlim([0, 1])
ax.set_ylim([0, 1])

ax.legend(loc='upper center')

ax.set_xlabel('Allele frequency, ' + r'$f$', fontsize=10)
ax.set_ylabel('Weighted function', fontsize=10)


fig.savefig(config.analysis_directory + "plot_weight.png", format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
