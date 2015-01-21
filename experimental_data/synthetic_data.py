from collections import defaultdict
import csv
from order.hist import Histogram


def read_synthetic_data_raw(filename='hist_by_ms_len_as_0_sum.tab'):
    hist_per_length_per_cycle = defaultdict(lambda: defaultdict(list))
    with open(filename, 'rb') as f:
        rdr = csv.DictReader(f, dialect='excel-tab')
        for row in rdr:
            hist_values = {key: int(row[str(key)]) for key in range(3, 100)}
            hist_per_length_per_cycle[int(row['len'])][int(row['cycles'])] = hist_values
        return hist_per_length_per_cycle


def to_hist_objects(hist_per_length_per_cycle, normalize=False, nsamples=None, trunc=False, cut_peak=False):
    hist_object_per_length_per_cycle = defaultdict(lambda: defaultdict(list))
    for ms_len in hist_per_length_per_cycle:
        for cycle in hist_per_length_per_cycle[ms_len]:
            hist_object = Histogram(hist_per_length_per_cycle[ms_len][cycle], normalize=normalize, nsamples=nsamples, trunc=trunc, cut_peak=cut_peak)
            hist_object_per_length_per_cycle[ms_len][cycle] = hist_object
    return hist_object_per_length_per_cycle


def transpose_hist_objects(hist_per_length_per_cycle):
    trans_hist_per_length_per_cycle = defaultdict(lambda: defaultdict(list))
    for ms_len in hist_per_length_per_cycle:
        for cycle in hist_per_length_per_cycle[ms_len]:
            trans_hist_per_length_per_cycle[ms_len][cycle] = hist_per_length_per_cycle[ms_len][cycle] - ms_len
    return trans_hist_per_length_per_cycle


def get_hists_pairs_by_len(hist_per_length_per_cycle, length=[5, 10, 15, 20, 25, 30]):
    hist_pairs = []
    for ms_len in length:
        syn_hist_26 = hist_per_length_per_cycle[ms_len][26]
        syn_hist_47 = hist_per_length_per_cycle[ms_len][47]
        hist_pairs.append((ms_len, syn_hist_26, syn_hist_47))
    return hist_pairs


def get_hists_pairs(filename='hist_by_ms_len_as_0_sum.tab', normalize=False, nsamples=None, trunc=False, cut_peak=False, length=[5, 10, 15, 20, 25, 30]):
    hist_per_length_per_cycle = read_synthetic_data_raw(filename=filename)
    hist_object_per_length_per_cycle = to_hist_objects(hist_per_length_per_cycle, normalize=normalize, nsamples=nsamples, trunc=trunc, cut_peak=cut_peak)
    return get_hists_pairs_by_len(hist_object_per_length_per_cycle, lengths=lengths)


def get_transposed_hists_pairs(filename='hist_by_ms_len_as_0_sum.tab', normalize=False, nsamples=None, trunc=False, cut_peak=False, length=[5, 10, 15, 20, 25, 30]):
    hist_per_length_per_cycle = read_synthetic_data_raw(filename=filename)
    hist_object_per_length_per_cycle = to_hist_objects(hist_per_length_per_cycle, normalize=normalize, nsamples=nsamples, trunc=trunc, cut_peak=cut_peak)
    trans_hist_per_length_per_cycle = transpose_hist_objects(hist_object_per_length_per_cycle)
    return get_hists_pairs_by_len(trans_hist_per_length_per_cycle, lengths=lengths)