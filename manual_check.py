from experimental_data.synthetic_data import read_synthetic_data, get_hists_pairs
from order.calibration.score import distance_from_model, distance_from_model_across_lengths


def test_score_all(up_x3, up_x2, up_x, up_c, dw_x3, dw_x2, dw_x, dw_c, method='bin', cycles_pair=(11,38)):
    up = lambda x: up_x3*x**3 + up_x2*x**2 + up_x*x + up_c
    dw = lambda x: dw_x3*x**3 + dw_x2*x**2 + dw_x*x + dw_c
    hist_pairs = get_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab')
    msmodel = (up, dw, method)
    return distance_from_model_across_lengths(msmodel, hist_pairs, cycles_pair)


def test_score_single(up_x3, up_x2, up_x, up_c, dw_x3, dw_x2, dw_x, dw_c, ms_len, method='bin', cycles_pair=(11,38)):
    up = lambda x: up_x3*x**3 + up_x2*x**2 + up_x*x + up_c
    dw = lambda x: dw_x3*x**3 + dw_x2*x**2 + dw_x*x + dw_c
    hist_per_length_per_cycle = read_synthetic_data(filename='experimental_data/hist_by_ms_len_as_0_sum.tab')
    hist_pairs = [(ms_len, hist_per_length_per_cycle[ms_len][26], hist_per_length_per_cycle[ms_len][47]),]
    msmodel = (up, dw, method)
    return distance_from_model_across_lengths(msmodel, hist_pairs, cycles_pair)
