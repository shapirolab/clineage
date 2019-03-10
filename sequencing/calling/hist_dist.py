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


def logprob(real_hist, sim_hist):
    np.dot(np.log(real_hist._vec), sim_hist._vec)
    return


def helper_dot(hv1, hv2):
    h1_n = np.linalg.norm(hv1, ord=2)
    h2_n = np.linalg.norm(hv2, ord=2)
    dot = np.dot(hv1,hv2) / (h1_n*h2_n)
    return float(dot)


def derived_proportions_dot(rht, sim_hist, eps=0.05):
    (rhmr, rhmc, rhv_n) = rht
    xt = np.matrix([sim_hist.vh1, sim_hist.vh2])
    x = xt.transpose()
    beta = np.linalg.inv(xt*x)*(xt*rhmc)
    if min(beta)<0 or sum(beta)==0:
        if beta[0] <= 0:
            return helper_dot(rhmr, sim_hist.vh2), 0.0
        return helper_dot(rhmr, sim_hist.vh1), 1.0
    if beta[0] / sum(beta) < eps:
        return helper_dot(rhmr, sim_hist.vh2), 0.0
    if beta[0] / sum(beta) > 1 - eps:
        return helper_dot(rhmr, sim_hist.vh1), 1.0
    p = beta[0]/sum(beta)
    yhat = x*beta
    yhat_n = np.linalg.norm(yhat, ord=2)
    conf = np.dot(rhmr, yhat)/(yhat_n*rhv_n)
    return float(conf), float(p)


def get_distance_function_by_name(func_name):
    if func_name == 'con':
        return pop_dist_corr_numpy
    if func_name == 'dot':
        return dot_product
    if func_name == 'dba':
        return derived_proportions_dot
    if func_name == 'sub':
        return substruction
    if func_name == 'ml':
        return pop_dist_corr_numpy
