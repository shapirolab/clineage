from ..preprocessing import generate_hist
from ..hist_dist import pop_dist
from scipy import optimize
import numpy


def distance_from_model(syn_len, syn_hist, cycles, msmodel, distance_measure='con', **kwargs):
    """

    :param syn_len:
    :param syn_hist:
    :param cycles:
    :param msmodel: (up, dw, method)
    :param distance_measure:
    :return:
    """
    up, dw, method = msmodel
    model_hist = generate_hist(syn_len, cycles, method, up=up, dw=dw, **kwargs)
    for h in [model_hist, syn_hist]:
        h.sq_normalize()
        if len(h.keys()) == 0:
            print h
            raise
    return pop_dist(syn_hist, model_hist, method=distance_measure)


def distance_from_model_across_lengths(msmodel, syn_hist_list, cycles_pair, distance_measure='con', **kwargs):
    """

    :param msmodel: (up, dw, method)
    :param syn_hist_list: [(ms_len, syn_hist_26, syn_hist_47), ...]
    :param cycles_pair: (11,38)
    :return:
    """
    cycles_26, cycles_47 = cycles_pair
    score = 0.0
    for syn_len, syn_hist_26, syn_hist_47 in syn_hist_list:
        score += distance_from_model(syn_len, syn_hist_26, cycles_26, msmodel, distance_measure=distance_measure, **kwargs)
        score += distance_from_model(syn_len, syn_hist_47, cycles_47, msmodel, distance_measure=distance_measure, **kwargs)
    # print msmodel, score
    return score


def get_x_from_list(l, step):
    assert len(l) % step == 0
    for i in xrange(0, len(l), step):
        res = []
        for j in xrange(step):
            res.append(l[i+j])
        yield tuple(res)


def optimize_across_lengths(input_tuple):
    alg, sim, optimizer_method, hist_pairs, cycles_tup, bounds, initial_guess, iterations, optimizer_options = input_tuple

    def nmes(x):
        assert len(x) % 2 == 0
        up_params = list(x)[:len(x)/4]
        dw_params = list(x)[len(x)/4:]
        # ups = [numpy.poly1d([a, b, c]) for a, b, c in get_x_from_list(up_params, 3)]
        # dws = [numpy.poly1d([a, b, c]) for a, b, c in get_x_from_list(dw_params, 3)]
        ups = [numpy.poly1d([a, b]) for a, b in get_x_from_list(up_params, 2)]
        dws = [numpy.poly1d([a, b]) for a, b in get_x_from_list(dw_params, 2)]
        msmodel = (ups, dws, sim)
        return distance_from_model_across_lengths(msmodel, hist_pairs, cycles_tup, distance_measure=alg)

    minimizer_kwargs = dict(method=optimizer_method, bounds=bounds, options=optimizer_options)
    res = optimize.basinhopping(nmes, initial_guess, minimizer_kwargs=minimizer_kwargs, niter=iterations)
    return alg, sim, optimizer_method, cycles_tup, [hp[0] for hp in hist_pairs], res