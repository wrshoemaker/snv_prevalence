








        smoothed_rsquared_numerators = numpy.array( smoothed_rsquared_numerators )
        smoothed_rsquared_denominators = numpy.array( smoothed_rsquared_denominators )
        smoothed_counts = numpy.array( smoothed_counts )

        all_smoothed_rsquared_numerators = numpy.array( all_smoothed_rsquared_numerators )
        all_smoothed_rsquared_denominators = numpy.array( all_smoothed_rsquared_denominators )
        all_smoothed_counts = numpy.array( all_smoothed_counts )

        print(len(distances), len(rsquared_numerators))

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

        print('1', len(distances), len(rsquared_numerators))


        all_distances = smoothed_distances
        #all_distances = dmins
        all_rsquareds = all_smoothed_rsquared_numerators/(all_smoothed_rsquared_denominators)
        all_ns = all_smoothed_counts

        all_distances = all_distances[all_ns>0]
        all_rsquareds = all_rsquareds[all_ns>0]
        all_ns = all_ns[all_ns>0]


        print('2', len(distances), len(rsquared_numerators))
