

def predict_prevalence(species_to_run='all'):


    subject_sample_map = parse_HMP_data.parse_subject_sample_map()

    if species_to_run == 'all':
        species_to_run = good_species_list

    else:
        species_to_run = [species_to_run]


    prevalence_dict_mapgd = {}
    #for species_name in good_species_list:
    for species_name in species_to_run:

        sys.stderr.write("%s\n" % species_name)

        gene_location_name_dict = get_gene_location_name_dict(species_name)

        #sys.stderr.write("Loading whitelisted genes...\n")
        #core_genes = core_gene_utils.parse_core_genes(species_name)

        # Holds panel wide prevalence for each species
        #os.system('mkdir -p %ssnp_prevalences' % config.data_directory)
        #snp_file_path = "%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name)

        #if os.path.isfile(snp_file_path) == False:
        #    continue

        parameter_dict_path = "%sparameter_dict_mapgd/%s.dat" % (config.data_directory, species_name)
        pi_dict_path = "%spi_annotated_snps_mapgd/%s.dat" % (config.data_directory, species_name)

        if os.path.isfile(parameter_dict_path) == False:
            continue

        if os.path.isfile(pi_dict_path) == False:
            continue

        sys.stderr.write("Loading pi dict...\n")
        pi_dict = load_pi_dict(species_name)
        sys.stderr.write("Done!\n")

        #sys.stderr.write("Loading parameter dict...\n")
        #parameter_dict = load_parameter_dict(species_name)
        #sys.stderr.write("Done!\n")

        #mapgd_samples_to_keep_idx = numpy.asarray([numpy.where(samples_all == sample)[0][0] for sample in mapgd_samples])
        things_to_measure = ['f_max_all', 'observed_prevalence_all', 'f_mean_all', 'predicted_prevalence_all', 'predicted_f_mean_all', 'f_max_mapgd', 'observed_prevalence_mapgd', 'predicted_prevalence_mapgd', 'f_mean_mapgd', 'predicted_f_mean_mapgd', 'f_no_zeros_mapgd', 'f_var_mapgd', 'predicted_f_var_mapgd', 'n_non_zero_frequencies', 'sites', 'mean_coverage_mapgd', 'predicted_prevalence_mapgd_slm', 'observed_prevalence_mapgd_slm']
        # observed_prevalence_mapgd_slm_cov_20,  predicted_prevalence_mapgd_slm_cov_20
        #'f_mean_no_zeros_mapgd', 'f_var_no_zeros_mapgd', 'n_hosts'
        prevalence_dict_mapgd[species_name] = {}
        for clade_type in clade_types:
            prevalence_dict_mapgd[species_name][clade_type] = {}
            for pi_type in ['pi_exclude_boundary', 'pi_include_boundary']:
                prevalence_dict_mapgd[species_name][clade_type][pi_type] = {}
                for allowed_variant_type in allowed_variant_types:
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][allowed_variant_type] = {}
                    for t in things_to_measure:
                        prevalence_dict_mapgd[species_name][clade_type][pi_type][allowed_variant_type][t] = []


        # now go back and calculate
        for clade_type in clade_types:

            sys.stderr.write("%s\n" % clade_type)

            sys.stderr.write("Loading parameter dict...\n")
            parameter_dict_clade = load_parameter_dict(species_name, clade_type)
            sys.stderr.write("Done!\n")
            # load dict_with_all_data
            # results from naive allele frequency estimates (NOT mapgd)
            #prevalence_dict_all = {}
            #intermediate_filename_template_prevalence = config.data_directory+"predicted_observed_prevalence/%s_%s_%s.dat"
            #for pi_type in ['pi_exclude_boundary', 'pi_include_boundary']:
            #    with open(intermediate_filename_template_prevalence % (species_name, clade_type, pi_type), 'rb') as handle:
            #        prevalence_all_dict = pickle.load(handle)
            #        prevalence_dict_all[pi_type] = prevalence_all_dict


            for key, value in parameter_dict_clade.items():

                # use the depth of coverage used by MAPGD
                depths_i = value['depths']
                f_i = value['f_mapgd']
                samples_i = value['samples']
                location_aa_i = value['location_aa']
                f_max_i = value['f_max']
                f_mean_i = value['f_mean']
                f_var_i = value['f_var']
                f_no_zeros_i = value['f_no_zeros_mapgd']
                f_mean_no_zeros_i = value['f_mean_no_zeros']
                f_var_no_zeros_i = value['f_var_no_zeros']

                for pi_type in ['pi_exclude_boundary', 'pi_include_boundary']:

                    #if key in prevalence_dict_all[pi_type][location_aa_i]['sites']:
                    #    prevalence_idx = prevalence_dict_all[pi_type][location_aa_i]['sites'].index(key)

                    #else:
                    #    continue

                    pi_filter = []
                    depths_i_filter = []
                    f_i_filter = []

                    for sample_idx, sample in enumerate(samples_i):
                        if sample in pi_dict[location_aa_i]:
                            if pi_type in pi_dict[location_aa_i][sample]:
                                #mapgd_samples_to_keep.append(sample)
                                depths_i_idx = depths_i[sample_idx]

                                #if depths_i_idx < D_min:
                                #    continue

                                pi_filter.append(pi_dict[location_aa_i][sample][pi_type])
                                #depths_i_filter.append(depths_i[sample_idx])
                                depths_i_filter.append(depths_i_idx)
                                f_i_filter.append(f_i[sample_idx])


                    pi_filter = numpy.asarray(pi_filter)
                    depths_i_filter = numpy.asarray(depths_i_filter)
                    f_i_filter = numpy.asarray(f_i_filter)

                    if len(pi_filter) < min_n_samples:
                        continue


                    #f_max_i_min = min([f_max_i, 1-f_max_i])

                    predicted_prevalence = 1 - numpy.mean( (1 +  f_max_i * depths_i_filter) ** (-1*pi_filter) )
                    predicted_f_mean = f_max_i * numpy.mean(pi_filter)
                    predicted_f_var = (f_max_i**2) * numpy.mean(pi_filter)

                    mean_coverage_mapgd = numpy.mean(depths_i_filter)
                    observed_prevalence = sum(f_i_filter>0)/len(f_i_filter)

                    predicted_prevalence_slm, observed_prevalence_slm = predict_prevalence_slm(f_i_filter, depths_i_filter)
                    # error SLM

                    error_slm = numpy.absolute(observed_prevalence_slm - predicted_prevalence_slm) / observed_prevalence_slm
                    #error_f = numpy.absolute(observed_prevalence - predicted_prevalence) / observed_prevalence


                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_max_mapgd'].append(f_max_i)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['observed_prevalence_mapgd'].append(value['observed_prevalence'])
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['predicted_prevalence_mapgd'].append(predicted_prevalence)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_mean_mapgd'].append(f_mean_i)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['predicted_f_mean_mapgd'].append(predicted_f_mean)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_var_mapgd'].append(f_var_i)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['predicted_f_var_mapgd'].append(predicted_f_var)

                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_no_zeros_mapgd'].extend(f_no_zeros_i)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['n_non_zero_frequencies'].append(len(f_no_zeros_i))
                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['n_hosts'].append(len(n_hosts))

                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_mean_no_zeros_mapgd'].append(f_mean_no_zeros_i)
                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_var_no_zeros_mapgd'].append(f_var_no_zeros_i)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['mean_coverage_mapgd'].append(mean_coverage_mapgd)

                    # predicted_prevalence_mapgd_slm
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['predicted_prevalence_mapgd_slm'].append(predicted_prevalence_slm)
                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['observed_prevalence_mapgd_slm'].append(observed_prevalence_slm)




                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_max_all'].append(prevalence_dict_all[pi_type][location_aa_i]['f_max'][prevalence_idx])
                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['observed_prevalence_all'].append(prevalence_dict_all[pi_type][location_aa_i]['observed_prevalence'][prevalence_idx])
                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['predicted_prevalence_all'].append(prevalence_dict_all[pi_type][location_aa_i]['predicted_prevalence'][prevalence_idx])
                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['f_mean_all'].append(prevalence_dict_all[pi_type][location_aa_i]['f_mean'][prevalence_idx])
                    #prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['predicted_f_mean_all'].append(prevalence_dict_all[pi_type][location_aa_i]['predicted_f_mean'][prevalence_idx])

                    prevalence_dict_mapgd[species_name][clade_type][pi_type][location_aa_i]['sites'].append(key)



    sys.stderr.write("Saving prevalence dict...\n")
    intermediate_filename = config.data_directory+"predicted_prevalence_mapgd_test.dat"
    with open(intermediate_filename, 'wb') as handle:
        pickle.dump(prevalence_dict_mapgd, handle, protocol=pickle.HIGHEST_PROTOCOL)
