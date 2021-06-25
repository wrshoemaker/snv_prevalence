import sys
import numpy

import config
import parse_midas_data
import diversity_utils
import sample_utils
import calculate_linkage_disequilibria

#import calculate_snv_distances
import figure_utils
from math import log10,ceil
#import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.random import randint, multinomial
import parse_HMP_data

#focal_speciess = ['Bacteroides_vulgatus_57955', 'Akkermansia_muciniphila_55290', 'Bacteroides_uniformis_57318']
focal_speciess = ['Bacteroides_vulgatus_57955']

species_name='Bacteroides_vulgatus_57955'

focal_colors = ['b','g']

variant_types = ['4D','1D']

color_dict = {'4D':'b', '1D':'r'}

lows = [0, 0.05, 0.1]#, 0.3, 0.5, 0.7]

#lines_styles = ['-', (0, (5, 1)), (0, (5, 10)), (0, (3, 5, 1, 5)), (0, (1, 1)),  (0, (1, 10)) ]
#lines_styles = ['-', '--', '-.', ':']
lines_styles = ['-', '--',  ':']

passed_species = []
sample_sizes = []


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize



if species_name in focal_speciess:
    focal_species_idx = focal_speciess.index(species_name)
    color = focal_colors[focal_species_idx]
    linewidth=1
    zorder=2
    alpha=1
else:
    color = 'r'
    alpha =0.3
    linewidth=0.5
    zorder = 1




sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

sys.stderr.write("Loading haploid samples...\n")
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

sys.stderr.write("Calculating unique samples...\n")

snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]

# Load inconsistency data
#sys.stderr.write("(core genes only...)\n")
#snv_distance_map = calculate_snv_distances.load_snv_distance_map(species_name)

#ds = numpy.logspace(log10(3e-05),log10(3e-02),50) # 15 points are plotted
#total_snps = numpy.zeros_like(ds)
#inconsistent_snps = numpy.zeros_like(ds)


def calculate_ld_map(species_name, low=float(0), high=float(1)):

    # Load precomputed LD
    ld_map = calculate_linkage_disequilibria.load_ld_map(species_name, low=low, high=high)
    ld_dict = {}
    #print(ld_map.keys())

    for variant_type in variant_types:

        ld_dict[variant_type] = {}


        if len(ld_map)>0:

            distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerators, control_rsquared_denominators, control_n, pi = ld_map[('largest_clade',variant_type)]

            if True:
                passed_species.append(species_name)
                sample_sizes.append(len(snp_samples))
            else:
                sys.stderr.write("%s intergene LD too high: %g (%g)\n" % (species_name, control_rsquared, rsquareds[0]))


        all_distances, all_rsquared_numerators, all_rsquared_denominators, all_ns, all_intergene_distances, all_intergene_rsquared_numerators, all_intergene_rsquared_denominators, all_intergene_ns, all_control_rsquared_numerator, all_control_rsquared_denominator, all_control_n, all_pi = ld_map[('all',variant_type)]
        all_control_rsquared = all_control_rsquared_numerator/all_control_rsquared_denominator

        distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerator, control_rsquared_denominator, control_n, pi = ld_map[('largest_clade',variant_type)]
        control_rsquared = control_rsquared_numerator/control_rsquared_denominator

        # smooth this stuff:
        smoothed_distances = distances
        window_width = 10**(0.1)

        dmins = smoothed_distances/(window_width**0.5)
        dmaxs = smoothed_distances*(window_width**0.5)


        smoothed_rsquared_numerators = []
        smoothed_rsquared_denominators = []
        smoothed_counts = []

        all_smoothed_rsquared_numerators = []
        all_smoothed_rsquared_denominators = []
        all_smoothed_counts = []

        print(len(zip(dmins,dmaxs)))

        for dmin,dmax in zip(dmins,dmaxs):

            binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
            binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
            binned_counts = ns[(distances>=dmin)*(distances<=dmax)]

            #print(binned_counts.sum())

            if binned_counts.sum()==0:
                continue

            smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
            smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
            smoothed_counts.append( binned_counts.sum() )

            binned_numerators = all_rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
            binned_denominators = all_rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
            binned_counts = all_ns[(distances>=dmin)*(distances<=dmax)]
            all_smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
            all_smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
            all_smoothed_counts.append( binned_counts.sum() )



        smoothed_rsquared_numerators = numpy.array( smoothed_rsquared_numerators )
        smoothed_rsquared_denominators = numpy.array( smoothed_rsquared_denominators )
        smoothed_counts = numpy.array( smoothed_counts )

        all_smoothed_rsquared_numerators = numpy.array( all_smoothed_rsquared_numerators )
        all_smoothed_rsquared_denominators = numpy.array( all_smoothed_rsquared_denominators )
        all_smoothed_counts = numpy.array( all_smoothed_counts )

        early_distances = distances[distances<101]
        early_rsquareds = rsquared_numerators[distances<101]*1.0/rsquared_denominators[distances<101]
        early_ns = ns[distances<101]

        early_distances = early_distances[early_ns>0.5]
        early_rsquareds = early_rsquareds[early_ns>0.5]
        early_ns = early_ns[early_ns>0.5]

        distances = smoothed_distances
        rsquareds = smoothed_rsquared_numerators/(smoothed_rsquared_denominators)
        ns = smoothed_counts

        distances = distances[ns>0]
        rsquareds = rsquareds[ns>0]
        ns = ns[ns>0]

        #print(len(rsquared_numerators), len(distances), len(ns))


        all_distances = smoothed_distances
        #all_distances = dmins
        all_rsquareds = all_smoothed_rsquared_numerators/(all_smoothed_rsquared_denominators)
        all_ns = all_smoothed_counts

        all_distances = all_distances[all_ns>0]
        all_rsquareds = all_rsquareds[all_ns>0]
        all_ns = all_ns[all_ns>0]


        if (species_name in focal_speciess):
            focal_species_idx = focal_speciess.index(species_name)
            #example_axis=focal_example_axes[focal_species_idx]
            example_idx = -1
            color=focal_colors[focal_species_idx]
        else:
            example_idx = supplemental_focal_species.index(species_name)
            #example_axis = example_axes[example_idx]
            color=focal_colors[0]

        num_bootstraps = 10

        bootstrapped_sigmasquareds = [] # will eventually be a matrix where first index is window_idx and second index is bootstrap index (sorted from lowest to highest)
        # Estimate bootstrap intervals for focal species only
        for dmin,dmax in zip(dmins,dmaxs):
            #print(len(rsquared_numerators), len(distances))
            binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
            binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
            binned_counts = ns[(distances>=dmin)*(distances<=dmax)]

            total_pairs = binned_counts.sum()

            upper_rsquareds = []
            lower_rsquareds = []

            if total_pairs>0:

                if True: #len(binned_counts)>1:
                    #print total_pairs
                    #print binned_counts
                    ps = binned_counts*1.0/total_pairs

                    window_bootstrapped_countss = multinomial(total_pairs,ps,size=num_bootstraps)

                    #print window_bootstrapped_countss.shape
                    window_bootstrapped_numerators = (window_bootstrapped_countss*binned_numerators[None,:]).sum(axis=1)*1.0/total_pairs
                    window_bootstrapped_denominators = (window_bootstrapped_countss*binned_denominators[None,:]).sum(axis=1)*1.0/total_pairs

                    window_bootstrapped_sigmasquareds = window_bootstrapped_numerators/window_bootstrapped_denominators

                    #print window_bootstrapped_sigmasquareds.shape
                    window_bootstrapped_sigmasquareds.sort()

                    bootstrapped_sigmasquareds.append(window_bootstrapped_sigmasquareds)

                    #print total_pairs

                else:
                    bootstrapped_sigmasquareds.append([binned_numerators/binned_denominators]*num_bootstraps)

            else:

                bootstrapped_sigmasquareds.append([-1]*num_bootstraps)


        upper_rsquareds = numpy.array([bootstrapped_sigmasquareds[window_idx][int(num_bootstraps*0.95)] for window_idx in range(0,len(bootstrapped_sigmasquareds))])
        lower_rsquareds = numpy.array([bootstrapped_sigmasquareds[window_idx][int(num_bootstraps*0.05)] for window_idx in range(0,len(bootstrapped_sigmasquareds))])

        #print upper_rsquareds-lower_rsquareds

        good_distances = (upper_rsquareds>=-0.5)*(lower_rsquareds>=-0.5)

        #theory_ls = numpy.logspace(0,log10(distances[-1]),100)
        #theory_NRs = theory_ls/200.0
        #theory_rsquareds = (10+2*theory_NRs)/(22+26*theory_NRs+4*theory_NRs*theory_NRs)
        ld_dict[variant_type]['all_distances'] = all_distances
        ld_dict[variant_type]['distances'] = distances
        ld_dict[variant_type]['good_distances'] = good_distances
        ld_dict[variant_type]['early_distances'] = early_distances

        ld_dict[variant_type]['all_control_rsquared'] = all_control_rsquared

        ld_dict[variant_type]['rsquareds'] = rsquareds
        ld_dict[variant_type]['all_rsquareds'] = all_rsquareds
        ld_dict[variant_type]['lower_rsquareds'] = lower_rsquareds
        ld_dict[variant_type]['upper_rsquareds'] = upper_rsquareds
        ld_dict[variant_type]['early_rsquareds'] = early_rsquareds
        ld_dict[variant_type]['control_rsquared'] = control_rsquared

    return ld_dict



#fig, ax = plt.subplots(figsize=(12,4))


fig, ax = plt.subplots(figsize=(12,4))
fig.subplots_adjust(bottom= 0.15)

for low_idx, low in enumerate(lows):

    ax = plt.subplot2grid((1, 3), (0, low_idx))

    ld_dict = calculate_ld_map(species_name, low=low, high=float(1))
    lines_style = lines_styles[low_idx]

    ax.set_title('Largest clade')



    for variant_type in variant_types:
        ld_dict_variant = ld_dict[variant_type]

        #print(ld_dict_variant['distances'])

        good_distances = ld_dict_variant['good_distances']


        #ax.loglog(ld_dict_variant['all_distances'], ld_dict_variant['all_rsquareds'],'-',color='0.7',label='All samples, %s' % variant_type)
        #ax.fill_between(numpy.array([3.3e03,1e04]),numpy.array([1e-02,1e-02]), numpy.array([1,1]),color='w',zorder=20)

        #ax.loglog([ld_dict_variant['all_distances'][-1],6e03], [ld_dict_variant['all_rsquareds'][-1], ld_dict_variant['all_control_rsquared']],':',color='0.7',zorder=21)
        #ax.loglog([6e03], [ld_dict_variant['all_control_rsquared']],'o',color='0.7',markersize=3,markeredgewidth=0,zorder=21)

        #print(('%f' % low).rstrip('0').rstrip('.') )

        ax.fill_between(ld_dict_variant['distances'][good_distances], ld_dict_variant['lower_rsquareds'][good_distances], ld_dict_variant['upper_rsquareds'][good_distances], linewidth=0, color=color,alpha=0.5)
        ax.loglog(ld_dict_variant['distances'], ld_dict_variant['rsquareds'], lines_styles[0], color=color_dict[variant_type],label='%s $f_{min}=$%s' % (variant_type, ('%f' % low).rstrip('0').rstrip('.') ), alpha=0.7 )
        ax.loglog(ld_dict_variant['early_distances'], ld_dict_variant['early_rsquareds'],'o',color=color_dict[variant_type],markersize=2,markeredgewidth=0,alpha=0.5)

        ax.loglog([ld_dict_variant['distances'][-1],6e03], [ld_dict_variant['rsquareds'][-1], ld_dict_variant['control_rsquared']],':',color=color_dict[variant_type],zorder=21)
        ax.loglog([6e03], [ld_dict_variant['control_rsquared']],'o',color=color_dict[variant_type],markersize=3,markeredgewidth=0,zorder=21)



        #ax.loglog(all_distances, all_rsquareds,'-',color='0.7',label='All samples')
        #ax.fill_between(numpy.array([3.3e03,1e04]),numpy.array([1e-02,1e-02]), numpy.array([1,1]),color='w',zorder=20)

        #ax.loglog([all_distances[-1],6e03], [all_rsquareds[-1], all_control_rsquared],':',color='0.7',zorder=21)
        #ax.loglog([6e03], [all_control_rsquared],'o',color='0.7',markersize=3,markeredgewidth=0,zorder=21)
        #ax.fill_between(distances[good_distances],lower_rsquareds[good_distances], upper_rsquareds[good_distances], linewidth=0, color=color,alpha=0.5)
        #ax.loglog(distances, rsquareds,'-',color=color,label='Largest clade')
        #ax.loglog(early_distances, early_rsquareds,'o',color=color,markersize=2,markeredgewidth=0,alpha=0.5)

        #ax.loglog([distances[-1],6e03], [rsquareds[-1], control_rsquared],':',color=color,zorder=21)
        #ax.loglog([6e03], [control_rsquared],'o',color=color,markersize=3,markeredgewidth=0,zorder=21)
        #ax.set_title( figure_utils.get_pretty_spe#ax.loglog(theory_ls, theory_rsquareds/theory_rsquareds[0]*3e-01,'k-',linewidth=0.3,zorder=0,label='Neutral')





        leg = ax.legend(loc='lower left',frameon=False, ) #title=figure_utils.get_pretty_species_name(species_name,include_number=False))
        leg._legend_box.align = "left"

        #example_axis.set_title(figure_utils.get_pretty_species_name(species_name, include_number=True),fontsize=5)

        line, = ax.loglog([ld_dict['4D']['distances'][-1],ld_dict['4D']['distances'][-1]],[1e-02,1],'k:')
        line.set_dashes((0.5,1))





        #example_idx
        #if example_idx>0:
        #    print("Setting ytick labels")
        #     ax.set_yticklabels([])

        ax.xaxis.get_major_ticks()[-2].label1.set_visible(False)
        ax.xaxis.get_major_ticks()[-2].tick1line.set_visible(False)

        for tick_idx in range(1,7):

            ax.xaxis.get_minor_ticks()[-tick_idx].tick1line.set_visible( False)
            ax.xaxis.get_minor_ticks()[-tick_idx].tick2line.set_visible( False)




        ax.set_ylabel('Linkage disequilibrium, $\sigma^2_d$' )


        ax.set_xlabel('Distance between SNVs, $\ell$')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_zorder(22)

        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        ax.set_xlim([2,1e04])
        ax.set_ylim([0.2e-01,1])

#ax.text(6e03,4e-03,'Genome-\nwide',     horizontalalignment='center',fontsize='5')

#ax.text(0.5, 1, 'Min freq. = %f, max freq. = %f' % (round(lower_threshold,3), round(upper_threshold,3)) , fontsize=14 )


fig.savefig("%s%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
