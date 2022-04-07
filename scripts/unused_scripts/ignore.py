

####### Old function
def calculate_pi_from_snps_fmax_cutoff(f_max_range):

    intermediate_filename_template = config.data_directory+"pi_annotated_snps_f_max_cutoff/%s.dat"

    for species_name in good_species_list:

        sys.stderr.write("Loading whitelisted genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)

        snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        if os.path.isfile(snp_file_path) == False:
            continue

        # Open post-processed MIDAS output
        snp_file =  bz2.BZ2File(snp_file_path, "r")
        line = snp_file.readline() # header
        items = line.split()[1:]
        samples = numpy.array([item.strip() for item in items])


        pi_dict = {}

        for f_max_range_i in f_max_range:

            pi_dict[f_max_range_i] = {}

            for allowed_variant_type in allowed_variant_types:
                pi_dict[f_max_range_i][allowed_variant_type] = {}

                for sample in samples:
                    pi_dict[f_max_range_i][allowed_variant_type][sample] = {}

                    pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_sum_include_boundary'] = 0
                    pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_include_boundary'] = 0

                    pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_sum_exclude_boundary'] = 0
                    pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_exclude_boundary'] = 0


        sys.stderr.write("Iterating though SNPs...\n")
        num_sites_processed = 0
        for line in snp_file:

            num_sites_processed+=1
            #print num_sites_processed
            if num_sites_processed%50000==0:
                sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))
                #if debug:
                #    break

            items = line.split()
            # Load information about site
            info_items = items[0].split("|")
            chromosome = info_items[0]
            location = long(info_items[1])
            gene_name = info_items[2]
            location_aa = info_items[3]

            # only examine core genes
            if gene_name not in core_genes:
                continue

            if location_aa not in allowed_variant_types:
                continue

            alts = []
            depths = []
            for item in items[1:]:
                subitems = item.split(",")
                alts.append(long(subitems[0]))
                depths.append(long(subitems[1]))
            alts = numpy.array(alts)
            depths = numpy.array(depths)
            refs = depths-alts


            # only try and predict prevalence if there's non-zero coverage for >= 20 hos
            if sum(depths>D_min) < min_n_samples:
                continue

            # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
            population_prevalence = ((alts>=refs)*(depths>D_min)).sum()
            population_freq = population_prevalence*1.0/(depths>D_min).sum()


            if population_freq>0.5:
                # alternate allele is in the majority
                # re-polarize for now
                alts,refs = refs,alts

            # Next calculate fraction of cohort where population minor allele is present at >=10% within-host frequency
            alt_threshold = numpy.ceil(depths*0.1)+0.5 #at least one read above 10%.

            snp_prevalence = ((alts>=alt_threshold)*(depths>D_min)).sum()
            snp_freq = snp_prevalence*1.0/(depths>D_min).sum()


            #if (population_freq==0) and (snp_freq==0):
            #    continue

            #if (population_freq==1) and (snp_freq==1):
            #    continue

            #if (snp_freq==0) or (snp_freq==1):
            #    continue


            depths_filter = depths[depths>D_min]

            freqs = alts[depths>D_min] / depths_filter
            samples_to_keep = samples[depths>D_min]


            f_max = max(freqs)
            # remove these sites because we don't examine them with the gamma
            if (f_max == 1) or (f_max == 0):
                continue

            for f_max_range_i in f_max_range:

                if f_max_range_i <= f_max:

                    for sample_i, freq_i in zip(samples_to_keep, freqs):

                        pi_dict[f_max_range_i][location_aa][sample_i]['pi_sum_include_boundary'] += 2*freq_i*(1-freq_i)
                        pi_dict[f_max_range_i][location_aa][sample_i]['n_sites_include_boundary'] += 1

                        if (freq_i != 1) and (freq_i != 0):

                            pi_dict[f_max_range_i][location_aa][sample_i]['pi_sum_exclude_boundary'] += 2*freq_i*(1-freq_i)
                            pi_dict[f_max_range_i][location_aa][sample_i]['n_sites_exclude_boundary'] += 1


        for f_max_range_i in f_max_range:
            for allowed_variant_type in allowed_variant_types:
                for sample in samples:

                    if pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_include_boundary'] > 0:
                        pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_include_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_include_boundary']

                    if pi_dict[f_max_range_i][allowed_variant_type][sample]['n_sites_exclude_boundary'] > 0:
                        pi_dict[f_max_range_i][allowed_variant_type][sample]['pi_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']


        intermediate_filename = intermediate_filename_template % species_name

        with open(intermediate_filename, 'wb') as handle:
            pickle.dump(pi_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
