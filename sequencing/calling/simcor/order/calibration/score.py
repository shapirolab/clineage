from ..hist_dist import pop_dist


def distance_from_model(syn_len, syn_hist, model_cycle, distance_measure):
    model_hist = model_cycle.get_hist_for_length(syn_len)
    return pop_dist(syn_hist, model_hist, method=distance_measure)


def distance_from_model_across_lengths(model, reference, distance_measure):
    """
    :param msmodel: (up, dw, method)
    :param syn_hist_list: [(ms_len, syn_hist_26, syn_hist_47), ...]
    :param cycles_pair: (11,38)
    :return:
    """
    syn_hist_list = reference.syn_hist_list
    cycles_dict = reference.cycles_dict

    score = 0.0
    for syn_len, syn_hist_0, syn_hist_26, syn_hist_47 in syn_hist_list:
        cycles_tup = cycles_dict[syn_len]
        if len(cycles_tup) == 2:
            cycles_26, cycles_47 = cycles_tup
            cycles_0 = 0
        elif len(cycles_tup) == 3:
            cycles_0, cycles_26, cycles_47 = cycles_tup
        else:
            raise ValueError('unexpected cycles tuple')
        model_0 = model.get_for_cycles(cycles_0)
        model_26 = model.get_for_cycles(cycles_26)
        model_47 = model.get_for_cycles(cycles_47)

        score += distance_from_model(syn_len, syn_hist_0, model_0,
            distance_measure=distance_measure)
        score += distance_from_model(syn_len, syn_hist_26, model_26,
            distance_measure=distance_measure)
        score += distance_from_model(syn_len, syn_hist_47, model_47,
            distance_measure=distance_measure)
    return score
