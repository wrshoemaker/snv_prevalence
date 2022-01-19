
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
