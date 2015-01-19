from ..preprocessing import generate_hist
from ..hist_dist import pop_dist


def distance_from_model(syn_len, syn_hist, cycles, msmodel, distance_measure='cor', **kwargs):
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
        h.normalize()
    return pop_dist(syn_hist, model_hist, method=distance_measure)


def distance_from_model_across_lengths(msmodel, syn_hist_list, cycles_pair):
    """

    :param msmodel: (up, dw, method)
    :param syn_hist_list: [(ms_len, syn_hist_26, syn_hist_47), ...]
    :param cycles_pair: (11,38)
    :return:
    """
    cycles_26, cycles_47 = cycles_pair
    score = 0.0
    for syn_len, syn_hist_26, syn_hist_47 in syn_hist_list:
        score += distance_from_model(syn_len, syn_hist_26, cycles_26, msmodel, distance_measure='cor')
        score += distance_from_model(syn_len, syn_hist_47, cycles_47, msmodel, distance_measure='cor')
    return score