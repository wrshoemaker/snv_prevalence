




def make_pi_and_parameter_dict(species_to_run='all'):

    subject_sample_map = parse_HMP_data.parse_subject_sample_map()

    if species_to_run == 'all':

        species_to_run = good_species_list

    else:

        species_to_run = [species_to_run]


    #for species_name in good_species_list:
    for species_name in species_to_run:

        sys.stderr.write("%s\n" % species_name)

        gene_location_name_dict = get_gene_location_name_dict(species_name)

        sys.stderr.write("Loading whitelisted genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)

        # Holds panel wide prevalence for each species
        #os.system('mkdir -p %ssnp_prevalences' % config.data_directory)
        snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        if os.path.isfile(snp_file_path) == False:
            continue

        # Open post-processed MIDAS output
        snp_file =  bz2.BZ2File(snp_file_path, "r")
        line = snp_file.readline() # header
        items = line.split()[1:]

        samples_all = numpy.array([item.strip() for item in items])

        sys.stderr.write("Loading MAPGD data...\n")
        mapgd_samples, mapgd_dict, chromosome_location_tuples = get_mapgd_data(species_name, samples_all, core_genes)
        mapgd_samples = numpy.asarray(mapgd_samples)
        mapgd_samples_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples])

        # get samples from largest clade
        snp_samples = diversity_utils.calculate_haploid_samples(species_name)
        if len(snp_samples) < min_sample_size:
            continue
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        # get clade idxs
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = numpy.array(dummy_samples)
        substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        coarse_grained_idxs, coarse_grained_cluster_list = clade_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold)
        coarse_grained_samples = snp_samples[coarse_grained_idxs]
        clade_sets = clade_utils.load_manual_clades(species_name)
        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)
        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])
        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]
        largest_clade_mapgd_samples = list(set(largest_clade_samples) & set(mapgd_samples))
        largest_clade_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in largest_clade_mapgd_samples])

        # no strains
        #mapgd_samples_no_strains = get_samples_no_strains(species_name, mapgd_samples)
        #mapgd_samples_no_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_no_strains])
        mapgd_samples_no_strains = get_samples_by_strain_status(species_name, mapgd_samples, strain_status=False)
        mapgd_samples_no_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_no_strains])


        # just strains
        #mapgd_samples_just_strains = get_samples_just_strains(species_name, mapgd_samples)
        #mapgd_samples_just_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_just_strains])
        mapgd_samples_just_strains = get_samples_by_strain_status(species_name, mapgd_samples, strain_status=True)
        mapgd_samples_just_strains_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples_just_strains])


        # get allelic states
        sys.stderr.write("Loading allelic states...\n")
        #allele_dict = get_allele_dict(species_name, mapgd_samples, chromosome_location_tuples)
        allele_dict = load_mapgd_alleles_dict(species_name)

        sys.stderr.write("Calculating SNP prevalences...\n")
        num_sites_processed = 0

        pi_dict = {}
        for allowed_variant_type in allowed_variant_types:
            pi_dict[allowed_variant_type] = {}
            for sample in mapgd_samples:
                pi_dict[allowed_variant_type][sample] = {}

                pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] = 0
                pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] = 0

                pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] = 0
                pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] = 0

        parameter_dict = {}
        for clade_type in clade_types:
            parameter_dict[clade_type] = {}

        # get chromosome and locations so you can go back and count alt alleles....
        for line in snp_file:

            num_sites_processed+=1
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

            if gene_name not in mapgd_dict:
                continue

            # Load alt and depth counts
            alts = []
            depths = []
            for item in items[1:]:
                subitems = item.split(",")
                alts.append(long(subitems[0]))
                depths.append(long(subitems[1]))
            alts = numpy.array(alts)
            depths = numpy.array(depths)
            refs = depths-alts

            # only try and predict prevalence if there's non-zero coverage for >= 20 hosts
            if sum(depths>D_min) < min_n_samples:
                continue

            # samples where mapgd was run and

            # only look at sites with mapgd data
            if location not in mapgd_dict[gene_name]:
                continue

            ref_allele = allele_dict[(chromosome, location)]['ref_allele']
            samples = mapgd_dict[gene_name][location]['samples']
            major_alleles = mapgd_dict[gene_name][location]['major_alleles']
            minor_alleles = mapgd_dict[gene_name][location]['minor_alleles']
            # frequencies are of MAJOR alleles
            frequencies = mapgd_dict[gene_name][location]['frequencies']

            #major_allele_counts = sorted(major_alleles, key=major_alleles.get, reverse=True)
            # choose the most common major allele as the reference
            frequencies_alts = []
            frequencies_alts_samples = []
            n_samples_alternate_major = 0
            for major_allele_idx, major_allele in enumerate(major_alleles):

                frequency = frequencies[major_allele_idx]
                # we want the frequencies of non reference alleles
                if major_allele == ref_allele:
                    frequency = 1 - frequency
                # count alternative alleles that are also major alles
                else:
                    n_samples_alternate_major += 1

                frequencies_alts.append(frequency)
                frequencies_alts_samples.append(samples[major_allele_idx])

                # add to pi, since 0 < frequency < 1
                pi_dict[location_aa][samples[major_allele_idx]]['pi_sum_include_boundary'] += 2*frequency*(1-frequency)
                pi_dict[location_aa][samples[major_allele_idx]]['pi_sum_exclude_boundary'] += 2*frequency*(1-frequency)
                pi_dict[location_aa][samples[major_allele_idx]]['n_sites_include_boundary'] += 1
                pi_dict[location_aa][samples[major_allele_idx]]['n_sites_exclude_boundary'] += 1


            for clade_type in clade_types:

                if clade_type == 'largest_clade':
                    alts_iter = alts[largest_clade_idx]
                    depths_iter = depths[largest_clade_idx]
                    refs_iter = refs[largest_clade_idx]
                    samples_to_save = samples_all[largest_clade_idx]
                    samples_clade_type = largest_clade_mapgd_samples


                elif clade_type == 'no_strains':
                    alts_iter = alts[mapgd_samples_no_strains_idx]
                    depths_iter = depths[mapgd_samples_no_strains_idx]
                    refs_iter = refs[mapgd_samples_no_strains_idx]
                    samples_to_save = samples_all[mapgd_samples_no_strains_idx]
                    samples_clade_type = mapgd_samples_no_strains


                elif clade_type == 'just_strains':
                    alts_iter = alts[mapgd_samples_just_strains_idx]
                    depths_iter = depths[mapgd_samples_just_strains_idx]
                    refs_iter = refs[mapgd_samples_just_strains_idx]
                    samples_to_save = samples_all[mapgd_samples_just_strains_idx]
                    samples_clade_type = mapgd_samples_just_strains


                else:
                    alts_iter = alts[mapgd_samples_idx]
                    depths_iter = depths[mapgd_samples_idx]
                    refs_iter = refs[mapgd_samples_idx]
                    samples_to_save = samples_all[mapgd_samples_idx]
                    samples_clade_type = mapgd_samples

                # get samples where MAPGD was run that are also in the clade type

                # original, overwrote by mistake
                #samples = [s for s in samples if s in samples_clade_type]
                samples_intersection = [s for s in samples if s in samples_clade_type]
                # get non zero frequencies for MAPGD samples that are in the clade type
                frequencies_alts_clade_type = []
                frequencies_alts_clade_type_samples = []
                for frequency_idx, frequency in enumerate(frequencies_alts):
                    if frequencies_alts_samples[frequency_idx] in samples_intersection:
                        frequencies_alts_clade_type.append(frequency)
                        frequencies_alts_clade_type_samples.append(frequencies_alts_samples[frequency_idx])

                # number of zeros
                #samples_with_absence = set(samples_clade_type) - set(samples_intersection)
                samples_with_absence = set(samples_clade_type) - set(frequencies_alts_clade_type_samples)
                #n_zeros_to_add = len(mapgd_samples) - len(samples)
                frequencies_alts_clade_type.extend([0]*len(samples_with_absence))
                frequencies_alts_clade_type = numpy.asarray(frequencies_alts_clade_type)
                # add samples back
                frequencies_alts_clade_type_samples.extend(samples_with_absence)

                # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
                population_prevalence = n_samples_alternate_major/len(samples_clade_type)
                reference = True
                if population_prevalence>0.5:
                    frequencies_alts_clade_type = 1 - frequencies_alts_clade_type
                    reference = False

                f_max = max(frequencies_alts_clade_type)

                #if f_max == float(1):
                #    continue

                f_mean = numpy.mean(frequencies_alts_clade_type)
                f_var = numpy.var(frequencies_alts_clade_type)

                prevalence = sum(frequencies_alts_clade_type>0)/len(frequencies_alts_clade_type)

                frequencies_alts_clade_type_no_zeros = frequencies_alts_clade_type[frequencies_alts_clade_type>0]

                if prevalence == 0:
                    continue

                if len(frequencies_alts_clade_type_no_zeros) > 3:
                    f_mean_no_zeros = numpy.mean(frequencies_alts_clade_type_no_zeros)
                    f_var_no_zeros = numpy.var(frequencies_alts_clade_type_no_zeros)
                else:
                    f_mean_no_zeros = float("nan")
                    f_var_no_zeros = float("nan")


                # go through all the samples where the allele was absent
                for s in samples_with_absence:
                    pi_dict[location_aa][s]['n_sites_include_boundary'] += 1

                parameter_dict[clade_type][(chromosome, location)] = {}

                parameter_dict[clade_type][(chromosome, location)]['f_max'] = f_max
                parameter_dict[clade_type][(chromosome, location)]['f_mean'] = f_mean
                parameter_dict[clade_type][(chromosome, location)]['f_var'] = f_var
                parameter_dict[clade_type][(chromosome, location)]['observed_prevalence'] = prevalence
                parameter_dict[clade_type][(chromosome, location)]['depths'] = depths_iter
                #parameter_dict[clade_type][(chromosome, location)]['samples'] = samples_to_save
                parameter_dict[clade_type][(chromosome, location)]['samples'] = frequencies_alts_clade_type_samples

                parameter_dict[clade_type][(chromosome, location)]['location_aa'] = location_aa
                parameter_dict[clade_type][(chromosome, location)]['f_no_zeros_mapgd'] = frequencies_alts_clade_type_no_zeros.tolist()

                parameter_dict[clade_type][(chromosome, location)]['f_mapgd'] = frequencies_alts_clade_type.tolist()

                parameter_dict[clade_type][(chromosome, location)]['f_mean_no_zeros'] = f_mean_no_zeros
                parameter_dict[clade_type][(chromosome, location)]['f_var_no_zeros'] = f_var_no_zeros


        snp_file.close()

        # calculate pi
        for allowed_variant_type in allowed_variant_types:
            for sample in mapgd_samples:
                if pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] > 0:
                    pi_dict[allowed_variant_type][sample]['pi_include_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_include_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_include_boundary']
                    pi_dict[allowed_variant_type][sample]['n_sites_include_boundary'] = pi_dict[allowed_variant_type][sample]['n_sites_include_boundary']

                if pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] > 0:
                    pi_dict[allowed_variant_type][sample]['pi_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['pi_sum_exclude_boundary'] / pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']
                    pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary'] = pi_dict[allowed_variant_type][sample]['n_sites_exclude_boundary']


        # save pi dict
        sys.stderr.write("Saving pi dict...\n")
        intermediate_filename_template_pi = config.data_directory+"pi_annotated_snps_mapgd/%s.dat"
        intermediate_filename_pi = intermediate_filename_template_pi % species_name

        with open(intermediate_filename_pi, 'wb') as handle:
            pickle.dump(pi_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


        # save parameter_dict
        sys.stderr.write("Saving parameter dict...\n")
        for clade_type in clade_types:
            intermediate_filename_template_parameter = config.data_directory+"parameter_dict_mapgd/%s_%s.dat"
            intermediate_filename_parameter = intermediate_filename_template_parameter % (species_name, clade_type)

            with open(intermediate_filename_parameter, 'wb') as handle:
                pickle.dump(parameter_dict[clade_type], handle, protocol=pickle.HIGHEST_PROTOCOL)
