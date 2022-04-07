import sys, os
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

focal_speciess = ['Bacteroides_vulgatus_57955', 'Akkermansia_muciniphila_55290']

species_name='Bacteroides_vulgatus_57955'

focal_colors = ['b','g']

variant_types = ['4D','1D']

color_dict = {'4D':'b', '1D':'r'}

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



color = 'blue'
linewidth=1
zorder=2
alpha=1

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species_list = parse_midas_data.parse_good_species_list()



def calculate_ld_dict(good_species_list, subject_sample_map, bin_width_exponent=0.1):

    ld_dict_all_species = {}

    for species_name in good_species_list:

        ld_directory = '%slinkage_disequilibria/' % (parse_midas_data.data_directory)
        intermediate_filename_template = '%s%s.txt.gz'
        intermediate_filename = intermediate_filename_template % (ld_directory, species_name)

        if os.path.exists(intermediate_filename) == False:
            continue

        sys.stderr.write("Loading haploid samples...\n")
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

        sys.stderr.write("Calculating unique samples...\n")

        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


        # Load precomputed LD
        ld_map = calculate_linkage_disequilibria.load_ld_map(species_name)
        ld_dict = {}

        for variant_type in variant_types:

            ld_dict[variant_type] = {}


            if len(ld_map)>0:

                distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerators, control_rsquared_denominators, control_n, pi = ld_map[('largest_clade',variant_type)]

                kegg_ids = ld_map[('largest_clade',variant_type)].keys()
                print(kegg_ids)

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
            #window_width = 10**(0.1)
            window_width = 10**(bin_width_exponent)

            dmins = smoothed_distances/(window_width**0.5)
            dmaxs = smoothed_distances*(window_width**0.5)

            smoothed_rsquared_numerators = []
            smoothed_rsquared_denominators = []
            smoothed_counts = []

            all_smoothed_rsquared_numerators = []
            all_smoothed_rsquared_denominators = []
            all_smoothed_counts = []

            for dmin,dmax in zip(dmins,dmaxs):
                binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
                binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
                binned_counts = ns[(distances>=dmin)*(distances<=dmax)]

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


            # adding this code below to fix error after filtering "distances"
            rsquared_numerators = rsquared_numerators[ns>0]
            rsquared_denominators = rsquared_denominators[ns>0]

            # then filter ns

            ns = ns[ns>0]

            all_distances = smoothed_distances
            #all_distances = dmins
            all_rsquareds = all_smoothed_rsquared_numerators/(all_smoothed_rsquared_denominators)
            all_ns = all_smoothed_counts

            all_distances = all_distances[all_ns>0]
            all_rsquareds = all_rsquareds[all_ns>0]
            all_ns = all_ns[all_ns>0]



            num_bootstraps = 10

            bootstrapped_sigmasquareds = [] # will eventually be a matrix where first index is window_idx and second index is bootstrap index (sorted from lowest to highest)
            # Estimate bootstrap intervals for focal species only
            for dmin,dmax in zip(dmins,dmaxs):
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


        ld_dict_all_species[species_name] = ld_dict


    return ld_dict_all_species




#good_species_list = ['Bacteroides_vulgatus_57955', 'Akkermansia_muciniphila_55290']

ld_dict_all_species = {}

for species_name in good_species_list:

    ld_directory = '%slinkage_disequilibria/' % (parse_midas_data.data_directory)
    intermediate_filename_template = '%s%s.txt.gz'
    intermediate_filename = intermediate_filename_template % (ld_directory, species_name)

    if os.path.exists(intermediate_filename) == False:
        continue

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

    print(species_name)


    # Load precomputed LD
    ld_map = calculate_linkage_disequilibria.load_ld_map(species_name)
    ld_dict = {}

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

        print(rsquared_numerators)
        print(distances)

        print(len(distances), len(rsquared_numerators))

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

        for dmin,dmax in zip(dmins,dmaxs):
            binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
            binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
            binned_counts = ns[(distances>=dmin)*(distances<=dmax)]
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


        all_distances = smoothed_distances
        #all_distances = dmins
        all_rsquareds = all_smoothed_rsquared_numerators/(all_smoothed_rsquared_denominators)
        all_ns = all_smoothed_counts

        all_distances = all_distances[all_ns>0]
        all_rsquareds = all_rsquareds[all_ns>0]
        all_ns = all_ns[all_ns>0]



        num_bootstraps = 10

        bootstrapped_sigmasquareds = [] # will eventually be a matrix where first index is window_idx and second index is bootstrap index (sorted from lowest to highest)
        # Estimate bootstrap intervals for focal species only
        for dmin,dmax in zip(dmins,dmaxs):
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


    #print(ld_dict['1D'])

    ld_dict_all_species[species_name] = ld_dict




fig, ax = plt.subplots(figsize=(4,4))

ax.axhline(y=1, ls=':', color='grey')


for species_name in ld_dict_all_species.keys():


    #print(len(['all_distances']))
    #print(len(ld_dict_all_species[species_name]['4D']['all_distances']))

    ld_dict = ld_dict_all_species[species_name]

    rsquareds_1D = ld_dict['1D']['rsquareds']
    rsquareds_4D = ld_dict['4D']['rsquareds']

    # remove zeros
    no_zeros_idx = rsquareds_4D > 0

    rsquareds_1D_no_zeros = ld_dict['1D']['rsquareds'][no_zeros_idx]
    rsquareds_4D_no_zeros = ld_dict['4D']['rsquareds'][no_zeros_idx]
    distances_1D_no_zeros = ld_dict['1D']['distances'][no_zeros_idx]
    distances_4D_no_zeros = ld_dict['4D']['distances'][no_zeros_idx]



    #good_distances = numpy.logical_and(ld_dict['1D']['good_distances'], ld_dict['4D']['good_distances'])
    #distances_good_distances = ld_dict['4D']['distances'][good_distances]


    #rsquareds_ratio = ld_dict['1D']['rsquareds'] / ld_dict['4D']['rsquareds']
    rsquareds_ratio_no_zeros = rsquareds_1D_no_zeros / rsquareds_4D_no_zeros

    rsquareds_ratio_no_zeros_cutoff = rsquareds_ratio_no_zeros[rsquareds_ratio_no_zeros>0.1]
    distances_4D_no_zeros_cutoff = distances_4D_no_zeros[rsquareds_ratio_no_zeros>0.1]

    #rsquareds_ratio_good_distances = rsquareds_ratio[good_distances]


    #rsquareds_ratio_good_distances_no_inf = rsquareds_ratio_good_distances[ (rsquareds_ratio_good_distances > -1* 1E308 ) & (rsquareds_ratio_good_distances < 1E308)]
    #distances_good_distances_no_inf = distances_good_distances[ (rsquareds_ratio_good_distances > -1* 1E308 ) & (rsquareds_ratio_good_distances < 1E308)]



    print(species_name, rsquareds_ratio_no_zeros_cutoff[rsquareds_ratio_no_zeros_cutoff<0.1])


    #rsquareds_ratio_filtered = rsquareds_ratio[ (rsquareds_ratio > -1* 1E308 ) & (rsquareds_ratio < 1E308)]


    #distances_filtered = ld_dict['4D']['distances'][ (rsquareds_ratio > -1* 1E308 ) & (rsquareds_ratio < 1E308)]

    #ax.loglog(ld_dict['4D']['distances'], rsquareds_ratio,'-',color='dodgerblue', alpha=0.7)
    ax.loglog(distances_4D_no_zeros_cutoff, rsquareds_ratio_no_zeros_cutoff,'-',color='dodgerblue', alpha=0.5)


    #ax.loglog(ld_dict_variant['all_distances'], ld_dict_variant['all_rsquareds'],'-',color='0.7',label='All samples, %s' % variant_type)

    rsquared_ratio_final = ld_dict['1D']['rsquareds'][-1] / ld_dict['4D']['rsquareds'][-1]
    control_rsquared_ratio = ld_dict['1D']['control_rsquared'] / ld_dict['4D']['control_rsquared']

    ax.loglog([ld_dict['4D']['distances'][-1],6e03], [rsquared_ratio_final, control_rsquared_ratio],':',color='dodgerblue',zorder=21)
    ax.loglog([6e03], [control_rsquared_ratio],'o',color='dodgerblue',markersize=3,markeredgewidth=0,zorder=21)




ax.set_ylabel('Linkage disequilibrium ratio, $\sigma^2_{d,N} / \sigma^2_{d,S}$' )

ax.set_xlabel('Distance between SNVs, $\ell$')
ax.text(6e03,2.6e-01,'Genome-\nwide',     horizontalalignment='center',fontsize='5')


from matplotlib.lines import Line2D

line = [Line2D([0], [0], color='dodgerblue', linestyle='-')]
label = ['Largest clade']
leg = ax.legend(line, label, loc='lower left',frameon=False, ) #title=figure_utils.get_pretty_species_name(species_name,include_number=False))
leg._legend_box.align = "left"

ax.set_title( "All species", fontsize=6)#,y=0.95)


fig.savefig("%s%s.png" % (config.analysis_directory, 'ld_ratio_all_species'), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



#for variant_type in variant_types:
#    print(variant_type)
#    ld_dict_variant = ld_dict[variant_type]


#    ax.loglog(ld_dict_variant['all_distances'], ld_dict_variant['all_rsquareds'],'-',color='0.7',label='All samples, %s' % variant_type)
#    ax.fill_between(numpy.array([3.3e03,1e04]),numpy.array([1e-02,1e-02]), numpy.array([1,1]),color='w',zorder=20)

#    ax.loglog([ld_dict_variant['all_distances'][-1],6e03], [ld_dict_variant['all_rsquareds'][-1], ld_dict_variant['all_control_rsquared']],':',color='0.7',zorder=21)
#    ax.loglog([6e03], [ld_dict_variant['all_control_rsquared']],'o',color='0.7',markersize=3,markeredgewidth=0,zorder=21)
#    ax.fill_between(ld_dict_variant['distances'][good_distances], ld_dict_variant['lower_rsquareds'][good_distances], ld_dict_variant['upper_rsquareds'][good_distances], linewidth=0, color=color,alpha=0.5)
#    ax.loglog(ld_dict_variant['distances'], ld_dict_variant['rsquareds'],'-',color=color_dict[variant_type],label='Largest clade, %s' % variant_type)
#    ax.loglog(ld_dict_variant['early_distances'], ld_dict_variant['early_rsquareds'],'o',color=color_dict[variant_type],markersize=2,markeredgewidth=0,alpha=0.5)

#    ax.loglog([ld_dict_variant['distances'][-1],6e03], [ld_dict_variant['rsquareds'][-1], ld_dict_variant['control_rsquared']],':',color=color_dict[variant_type],zorder=21)
#    ax.loglog([6e03], [ld_dict_variant['control_rsquared']],'o',color=color_dict[variant_type],markersize=3,markeredgewidth=0,zorder=21)







#line, = ax.loglog([ld_dict['4D']['distances'][-1],ld_dict['4D']['distances'][-1]],[1e-02,1],'k:')
#line.set_dashes((0.5,1))




#ax.xaxis.get_major_ticks()[-2].label1.set_visible(False)
#ax.xaxis.get_major_ticks()[-2].tick1line.set_visible(False)

#for tick_idx in range(1,7):

#    ax.xaxis.get_minor_ticks()[-tick_idx].tick1line.set_visible( False)
#    ax.xaxis.get_minor_ticks()[-tick_idx].tick2line.set_visible( False)




#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_zorder(22)

#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()

#ax.set_xlim([2,1e04])
#ax.set_ylim([1e-02,1])
