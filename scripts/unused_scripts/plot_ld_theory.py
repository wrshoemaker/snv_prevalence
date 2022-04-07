from __future__ import division
import  matplotlib.pyplot as plt

import numpy
import parse_midas_data

s = -0.0001
N = 10**6
mu = 0.00000001

#r = 0.00001

def ld_ratio(N, s, mu, r):

    return (1 - s + ((2*s*mu)/ ((2*mu)-s))) * (1 + ((240*s*N  + 36*s) / (45 - 240*s*N + 36*s + 138*N*(4*mu +r) + 128*((N*r)**2 ) )) )




r_list = [0.0001, 0.001, 0.01]

s_list = numpy.logspace(-4, -8, num=100)

ls_list = ['-', '--', ':']

fig, ax = plt.subplots(figsize=(4,4))

for r_i_idx, r_i in enumerate(r_list):

    ld_ratio_list = []

    for s_i in s_list:

        ld_ratio_list.append(ld_ratio(N, s_i, mu, r_i))


    ax.plot(s_list, ld_ratio_list, ls= ls_list[r_i_idx], label = r'$NR=%.2f$' % (N*r_i) )



ax.set_xlabel('Scaled selection coefficient, ' + r'$\left | Ns \right |$', fontsize=12)
ax.set_ylabel('Ratio of nonsynonymous to synonymous LD, ' + r'$\frac{\sigma^{2}_{d, N}}{\sigma^{2}_{d, S}}$', fontsize=10)

ax.legend(loc='upper right',  prop={'size': 8})

ax.set_xscale('log', basex=10)
#ax.set_yscale('log', basey=10)


fig.subplots_adjust(hspace=1.8, wspace=0.5) #hspace=0.3, wspace=0.5
fig_name = "%sld_ratio_theory.png" % parse_midas_data.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
