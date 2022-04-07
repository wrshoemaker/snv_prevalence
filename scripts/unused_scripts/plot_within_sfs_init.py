import matplotlib.pyplot as plt
import scipy.stats as stats

fig, ax = plt.subplots(figsize=(4,4))


for variant_type in variant_types:

    site_samples_matrix, sites, snp_samples, gene_site_map = calculate_sfs_labelled_matrix(species_name, allowed_variant_types=set([variant_type]), desired_samples=snp_samples, allowed_genes=allowed_genes)

    means = []
    variances = []

    for i in range(site_samples_matrix.shape[0]):
        frequencies = site_samples_matrix[i,:][site_samples_matrix[i,:] > 0]

        if len(frequencies) > 4:
            #print(sites[i] , frequencies)
            means.append(numpy.mean(frequencies))
            variances.append(numpy.var(frequencies))

    means_variances_zipped = [x for x in zip(means, variances) if x[0] < 0.2]
    means_cutoff = [x[0] for x in means_variances_zipped]
    variances_cutoff = [x[1] for x in means_variances_zipped]

    means = numpy.asarray(means)
    variances = numpy.asarray(variances)


    slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(means_cutoff), numpy.log10(variances_cutoff))

    ax.scatter(means, variances, alpha=0.05, c=variant_color[variant_type], label=variant_type)#, c='#87CEEB')


    x_log10_range =  numpy.linspace(min(numpy.log10(means_cutoff)) , max(numpy.log10(means_cutoff)) , 10000)
    y_log10_fit_range = 10 ** (slope*x_log10_range + intercept)
    #y_log10_null_range = 10 ** (slope_null*x_log10_range + intercept)


    t_value = (slope - slope_null)/std_err
    p_value = stats.t.sf(numpy.abs(t_value), len(means)-2)

    sys.stdout.write("Slope = %g, t = %g, P= %g\n" % (slope, t_value, p_value))

    ax.plot(10**x_log10_range, y_log10_fit_range,  c=variant_color[variant_type], lw=2.5, linestyle='-', zorder=2, label="OLS fit " + variant_type)
    #ax.plot(10**x_log10_range, y_log10_null_range, c='k', lw=2.5, linestyle='--', zorder=2, label="Taylor's law")

    if variant_type == '1D':
        y_ax = 0.9
    else:
        y_ax = 0.7

    ax.text(0.2, y_ax, r'$y \sim x^{{{}}}$'.format(str( round(slope, 3) )), fontsize=11, color=variant_color[variant_type], ha='center', va='center', transform=ax.transAxes  )


ax.set_xlim([0.05, 0.8])
ax.set_ylim([0.00001, 0.8])

ax.axvline(0.5, ls='--', lw=2, c='k')

ax.axvline(0.1, ls='--', lw=2, c='k')

ax.axvline(0.2, ls='--', lw=2, c='k')


ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)

ax.set_xlabel('Average SNP MAF', fontsize=12)
ax.set_ylabel('Variance of SNP MAFs', fontsize=10)


ax.legend(loc="lower right", fontsize=8)



fig.subplots_adjust(wspace=0.3, hspace=0.3)
fig.savefig(config.analysis_directory + "taylors_law.png", format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()


# now plot a


FFD_nonsynonymous = site_samples_matrix_nonsynonymous.flatten()
FFD_nonsynonymous = FFD_nonsynonymous[FFD_nonsynonymous>0]

FFD_synonymous = site_samples_matrix_synonymous.flatten()
FFD_synonymous = FFD_synonymous[FFD_synonymous>0]

fig, ax = plt.subplots(figsize=(4,4))
ax.hist(FFD_nonsynonymous, bins=40, weights=numpy.zeros_like(FFD_nonsynonymous) + 1. / len(FFD_nonsynonymous), color='r', label='1D', histtype='step')
ax.hist(FFD_synonymous, bins=40, weights=numpy.zeros_like(FFD_synonymous) + 1. / len(FFD_synonymous), color='g', label='4D', histtype='step')
ax.legend(loc="upper right", fontsize=8)


ax.set_xlabel('Minor alelle frequency', fontsize=12)
ax.set_ylabel('Fraction of SNPs', fontsize=10)
#fig.subplots_adjust(wspace=0.3, hspace=0.3)
fig.savefig(config.analysis_directory + "FFD.png", format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
