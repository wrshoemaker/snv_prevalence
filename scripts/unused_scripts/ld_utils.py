from __future__ import division

import sys, os
import numpy


# color gradient for taxa
def cmdscale(D):
    """
    Classical multidimensional scaling (MDS)

    Parameters
    ----------
    D : (n, n) array
        Symmetric distance matrix.

    Returns
    -------
    Y : (n, p) array
        Configuration matrix. Each column represents a dimension. Only the
        p dimensions corresponding to positive eigenvalues of B are returned.
        Note that each dimension is only determined up to an overall sign,
        corresponding to a reflection.

    e : (n,) array
        Eigenvalues of B.

    """
    # Number of points
    n = len(D)

    # Centering matrix
    H = numpy.eye(n) - numpy.ones((n, n))/n

    # YY^T
    B = -H.dot(D**2).dot(H)/2

    # Diagonalize
    evals, evecs = numpy.linalg.eigh(B)

    # Sort by eigenvalue in descending order
    idx   = numpy.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]

    # Compute the coordinates using positive-eigenvalued components only
    w, = numpy.where(evals > 0)
    L  = numpy.diag(numpy.sqrt(evals[w]))
    V  = evecs[:,w]
    Y  = V.dot(L)

    return Y, evals




def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in xrange(0, len(lst), n):
        yield lst[i:i + n]



def calculate_ld_4D_1D_ratio(ld_dict):

    ld_1D_distances = ld_dict['1D']['distances']
    ld_4D_distances = ld_dict['4D']['distances']

    #ld_1D_distances = ld_1D_distances[ld_dict['1D']['good_distances']]
    #ld_4D_distances = ld_4D_distances[ld_dict['4D']['good_distances']]

    # keep the distances for 1D and 4D


    distance_bins = numpy.linspace(3, 3000, num=1000, endpoint=True, retstep=3)[0]
    distance_bins_list = distance_bins.tolist()


    good_distances_1D = ld_dict['1D']['good_distances']
    good_distances_1D_list = good_distances_1D.tolist()

    good_distances_4D = ld_dict['4D']['good_distances']
    good_distances_4D_list = good_distances_4D.tolist()


    ld_1D_rsquareds = ld_dict['1D']['rsquareds']
    ld_4D_rsquareds = ld_dict['4D']['rsquareds']

    ld_1D_all_rsquareds = ld_dict['1D']['all_rsquareds']
    ld_4D_all_rsquareds = ld_dict['4D']['all_rsquareds']



    # find items missing from distances, hacky but it works
    if len(ld_1D_distances) != len(good_distances_1D):
        # identify missing bins
        difference_1D = numpy.setdiff1d(distance_bins,  ld_1D_distances)
        difference_1D_idx_list = [distance_bins_list.index(x) for x in difference_1D]
        difference_1D_idx = numpy.asarray(difference_1D_idx_list)
        good_distances_1D = numpy.delete(good_distances_1D, difference_1D_idx)


    # same for 4D....
    if len(ld_4D_distances) != len(good_distances_4D):
        # identify missing bins
        difference_4D = numpy.setdiff1d(distance_bins,  ld_4D_distances)
        difference_4D_idx_list = [distance_bins_list.index(x) for x in difference_4D]
        difference_4D_idx = numpy.asarray(difference_4D_idx_list)
        good_distances_4D = numpy.delete(good_distances_4D, difference_4D_idx)


    ld_1D_distances_list = ld_1D_distances.tolist()
    ld_4D_distances_list = ld_4D_distances.tolist()

    ld_filter_distances = numpy.intersect1d(ld_1D_distances[good_distances_1D], ld_4D_distances[good_distances_4D])
    # get distances to exclude, bullshit error keeps popping up
    ld_filter_distances = [x for x in ld_filter_distances if (ld_1D_distances_list.index(x) < len(ld_1D_rsquareds)) and (ld_1D_distances_list.index(x) < len(ld_1D_all_rsquareds)) and (ld_4D_distances_list.index(x) < len(ld_4D_rsquareds)) and (ld_4D_distances_list.index(x) < len(ld_4D_all_rsquareds))]
    ld_filter_distances = numpy.asarray(ld_filter_distances)



    ld_filter_distances_list = ld_filter_distances.tolist()

    ld_1D_filter_idx_list = [ld_1D_distances_list.index(x) for x in ld_1D_distances_list if x in ld_filter_distances_list]
    ld_4D_filter_idx_list = [ld_4D_distances_list.index(x) for x in ld_4D_distances_list if x in ld_filter_distances_list]

    ld_1D_filter_idx = numpy.asarray(ld_1D_filter_idx_list)
    ld_4D_filter_idx = numpy.asarray(ld_4D_filter_idx_list)

    # remove indexes greater than length of aray


    #ld_1D_filter_idx = numpy.where(ld_1D_distances == ld_filter_distances)[0]
    #ld_4D_filter_idx = numpy.where(ld_4D_distances == ld_filter_distances)[0]


    if (len(ld_1D_filter_idx) == 0) or (len(ld_4D_filter_idx) == 0):
        return [],[],[],0

    else:



        #max_position = min([len(ld_1D_rsquareds), len(ld_4D_rsquareds), len(ld_1D_all_rsquareds), len(ld_4D_all_rsquareds)])

        #ld_1D_filter_idx = ld_1D_filter_idx[ld_1D_filter_idx < max_position]
        #ld_4D_filter_idx = ld_4D_filter_idx[ld_4D_filter_idx < max_position]

        #ld_4D_filter_idx = ld_4D_filter_idx[ (ld_4D_filter_idx < len(ld_4D_rsquareds)) & (ld_4D_filter_idx < len(ld_4D_all_rsquareds))]

        #print(len(ld_1D_rsquareds), len(ld_1D_all_rsquareds), len(ld_4D_rsquareds), len(ld_4D_all_rsquareds))

        ld_1D_filter_rsquareds = ld_1D_rsquareds[ld_1D_filter_idx]
        ld_4D_filter_rsquareds = ld_4D_rsquareds[ld_4D_filter_idx]



        #print(len(ld_dict['1D']['all_rsquareds']), len(ld_dict['4D']['all_rsquareds']))
        #print(len(ld_1D_filter_idx), len(ld_4D_filter_idx), min(ld_1D_filter_idx), max(ld_1D_filter_idx), min(ld_4D_filter_idx), max(ld_4D_filter_idx))

        ld_1D_filter_all_rsquareds = ld_1D_all_rsquareds[ld_1D_filter_idx]
        ld_4D_filter_all_rsquareds = ld_4D_all_rsquareds[ld_4D_filter_idx]

        ld_1D_4D_ratio_filter_rsquareds = ld_1D_filter_rsquareds/ld_4D_filter_rsquareds

        ld_1D_4D_ratio_filter_all_rsquareds = ld_1D_filter_all_rsquareds/ld_4D_filter_all_rsquareds

        #rsquared_ratio_final = ld_1D_4D_ratio_filter_rsquareds[-1]
        control_rsquared_ratio = ld_dict['1D']['control_rsquared'] / ld_dict['4D']['control_rsquared']

        # get rid of the
        nan_idxs = numpy.isnan(ld_1D_4D_ratio_filter_rsquareds) + numpy.isnan(ld_1D_4D_ratio_filter_all_rsquareds)

        to_keep_idxs = numpy.logical_not(nan_idxs)

        ld_filter_distances = ld_filter_distances[to_keep_idxs]
        ld_1D_4D_ratio_filter_rsquareds = ld_1D_4D_ratio_filter_rsquareds[to_keep_idxs]
        ld_1D_4D_ratio_filter_all_rsquareds = ld_1D_4D_ratio_filter_all_rsquareds[to_keep_idxs]


        return ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds, ld_1D_4D_ratio_filter_all_rsquareds, control_rsquared_ratio


    #ax.loglog(ld_filter_distances, ld_1D_4D_ratio_filter_rsquareds,'-',color='b', alpha=0.5)
    #ax.loglog([  ld_filter_distances[-1],  6e03], [rsquared_ratio_final, control_rsquared_ratio],':',color='dodgerblue',zorder=21)
    #ax.loglog([6e03], [control_rsquared_ratio],'o',color='b',markersize=3,markeredgewidth=0,zorder=21)
