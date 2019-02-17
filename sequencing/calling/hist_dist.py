import numpy as np
from sequencing.calling.hist import get_lims
import sys


def maximum_likelihood(hist1, hist2):
    li, ri = get_lims(hist1, hist2)
    return -sum(np.log(hist1[x] if hist1[x] > 0 else sys.float_info.epsilon)*hist2[x] for x in range(li, ri))/sum(hist2._hist.values())


def pop_dist_corr_numpy(hist1, hist2):
    """
    1-correlation
    """
    li, ri = get_lims(hist1, hist2)
    pmat = np.corrcoef([hist1[x] for x in range(li, ri)], [hist2[x] for x in range(li, ri)])
    return 1-pmat[0][1]


def dot_product(hist1, hist2):
    li, ri = get_lims(hist1, hist2)
    h1 = hist1.copy()
    h2 = hist2.copy()
    h1.sq_normalize()
    h2.sq_normalize()
    return 1-np.dot([h1[x] for x in range(li, ri)], [h2[x] for x in range(li, ri)])


def dotv(hist1, hist2):
    h1_n = np.linalg.norm(hist1._vec, ord=2)
    h2_n = np.linalg.norm(hist2._vec, ord=2)
    dot = np.dot(hist1._vec,hist2._vec) / (h1_n*h2_n)
    return 1-dot


def substruction(hist1, hist2):
    h1_n = np.linalg.norm(hist1._vec, ord=1)
    h2_n = np.linalg.norm(hist2._vec, ord=1)
    return np.linalg.norm(hist1._vec-hist2._vec, ord=1)/(h1_n*h2_n)


def get_distance_function_by_name(func_name):
    if func_name == 'con':
        return pop_dist_corr_numpy
    if func_name == 'dot':
        return dot_product
    if func_name == 'sub':
        return substruction
    if func_name == 'ml':
        return pop_dist_corr_numpy
