from experimental_data.synthetic_data import read_transposed_synthetic_data, get_hists_pairs_by_len
from order.calibration.score import distance_from_model_across_lengths
from frogress import bar

def test_score_all(up_x3, up_x2, up_x, up_c, dw_x3, dw_x2, dw_x, dw_c, method='bin', cycles_pair=(11,38), distance_measure='cor', **kwargs):
    up = lambda x: up_x3*x**3 + up_x2*x**2 + up_x*x + up_c
    dw = lambda x: dw_x3*x**3 + dw_x2*x**2 + dw_x*x + dw_c
    hist_pairs = get_hists_pairs_by_len(read_transposed_synthetic_data(filename='experimental_data/hist_by_ms_len_as_0_sum.tab'))
    msmodel = (up, dw, method)
    return distance_from_model_across_lengths(msmodel, hist_pairs, cycles_pair, distance_measure=distance_measure, **kwargs)


def test_score_single(up_x3, up_x2, up_x, up_c, dw_x3, dw_x2, dw_x, dw_c, ms_len, method='bin', cycles_pair=(11,38), distance_measure='cor', **kwargs):
    up = lambda x: up_x3*x**3 + up_x2*x**2 + up_x*x + up_c
    dw = lambda x: dw_x3*x**3 + dw_x2*x**2 + dw_x*x + dw_c
    hist_per_length_per_cycle = read_transposed_synthetic_data(filename='experimental_data/hist_by_ms_len_as_0_sum.tab')
    hist_pairs = [(ms_len, hist_per_length_per_cycle[ms_len][26], hist_per_length_per_cycle[ms_len][47]),]
    msmodel = (up, dw, method)
    return distance_from_model_across_lengths(msmodel, hist_pairs, cycles_pair, distance_measure=distance_measure, **kwargs)


# In [3]: test_score_all(0, 0.00005,  -0.0009, 0.0036, 0, 0.00009, -0.00003, -0.0013, method='bin', cycles_pair=(11,38), distance_measure='cor')
# Out[3]: 0.12080802360769172

# In [4]: test_score_all(0, 0.00005,  -0.0009, 0.0036, 0, 0.00009, -0.00003, -0.0013, method='mat', cycles_pair=(11,38), distance_measure='cor')
# Out[4]: 0.13542764059100665

# In [5]: test_score_all(0, 0.00005,  -0.0009, 0.0036, 0, 0.00009, -0.00003, -0.0013, method='dyn', cycles_pair=(11*2,38*2), distance_measure='cor')
# Out[5]: 0.13206008114145618

def moddiff(d='con'):
    diffs = []
    for i in bar(range(100)):
        diffs.append(test_score_single(0, 0.00005,  -0.0009, 0.0036, 0, 0.00009, -0.00003, -0.0013, 15, method='mat', cycles_pair=(11,38), distance_measure=d, sample_depth=1000000) - test_score_single(0, 0.00005,  -0.0009, 0.0036, 0, 0.00009, -0.00003, -0.0013, 15, method='bin', cycles_pair=(11,38), distance_measure=d, sample_depth=1000000))
    return diffs

def count_winners(d='con'):
    diffs = moddiff(d='con')
    return sum([i>0 for i in diffs])/float(len(diffs))

for i in ['sub', 'sp', 'emd', 'cor', 'con', 'chi', 'pr', 'kl', 'ks']:
    print i, count_winners(d=i)