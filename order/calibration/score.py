from ..preprocessing import generate_mat_hist
from ..hist_dist import pop_dist
from order.optimize_probs import build_markovian_matrix
import numpy
from numpy.linalg import matrix_power


class hashable_poly1d(numpy.poly1d):
    def __hash__(self):
        return hash((tuple(self.coeffs), self.order, self.variable))


def distance_from_model(syn_len, syn_hist, model_mat, distance_measure='con', **kwargs):
    model_hist = generate_mat_hist(syn_len, model_mat, **kwargs)
    return pop_dist(syn_hist, model_hist, method=distance_measure)


def distance_from_model_across_lengths(steps, reference, distance_measure='con', **kwargs):
    """
    :param msmodel: (up, dw, method)
    :param syn_hist_list: [(ms_len, syn_hist_26, syn_hist_47), ...]
    :param cycles_pair: (11,38)
    :return:
    """
    syn_hist_list = reference.syn_hist_list
    cycles_tup = reference.cycles_tup
    if len(cycles_tup) == 2:
        cycles_26, cycles_47 = cycles_tup
        cycles_0 = 0
    elif len(cycles_tup) == 3:
        cycles_0, cycles_26, cycles_47 = cycles_tup
    else:
        raise ValueError('unexpected cycles tuple')
    base_mat = build_markovian_matrix(steps, **kwargs)
    mat_0 = matrix_power(base_mat, cycles_0)
    mat_26 = matrix_power(base_mat, cycles_26)
    mat_47 = matrix_power(base_mat, cycles_47)
    score = 0.0
    for syn_len, syn_hist_0, syn_hist_26, syn_hist_47 in syn_hist_list:
        score += distance_from_model(syn_len, syn_hist_0, mat_0, distance_measure=distance_measure, **kwargs)
        score += distance_from_model(syn_len, syn_hist_26, mat_26, distance_measure=distance_measure, **kwargs)
        score += distance_from_model(syn_len, syn_hist_47, mat_47, distance_measure=distance_measure, **kwargs)
    return score


def get_x_from_list(l, step):
    assert len(l) % step == 0
    for i in range(0, len(l), step):
        res = []
        for j in range(step):
            res.append(l[i+j])
        yield tuple(res)


def param_list_to_polynomes_mutation(x):
    assert len(x) % 2 == 0
    # up_params = x[:len(x)//2]
    # dw_params = x[len(x)//2:]
    # ups = [hashable_poly1d([a, b, c]) for a, b, c in get_x_from_list(up_params, 3)]
    # dws = [hashable_poly1d([a, b, c]) for a, b, c in get_x_from_list(dw_params, 3)]
    up_params = x[:len(x)//4]
    dw_params = x[len(x)//4:]
    ups = [hashable_poly1d([a, b]) for a, b in get_x_from_list(up_params, 2)]
    dws = [hashable_poly1d([a, b]) for a, b in get_x_from_list(dw_params, 2)]
    # ups = [Hashable_exp(a, b) for a, b in get_x_from_list(up_params, 2)]
    # dws = [Hashable_exp(a, b) for a, b in get_x_from_list(dw_params, 2)]
    
    def fs(n):
        return 1-sum([fu(n) for fu in ups])-sum([fd(n) for fd in dws])
    
    d = {0: fs}
    d.update({i+1: fu for i, fu in enumerate(ups)})
    d.update({-i-1: fd for i, fd in enumerate(dws)})
    
    return d


class Hashable_exp(object):
    def __init__(self, a, b):
        self._a = a
        self._b = b

    def __call__(self, x):
        return self._a*(numpy.exp(self._b*x)-1)
