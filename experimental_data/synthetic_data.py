from collections import defaultdict
import csv


def read_synthetic_data(filename='hist_by_ms_len_as_0_sum.tab'):
    hist_per_length_per_cycle = defaultdict(lambda: defaultdict(list))
    with open(filename, 'rb') as f:
        rdr = csv.DictReader(f, dialect='excel-tab')
        for row in rdr:
            hist_values = {key: int(row[str(key)]) for key in range(3, 100)}
            hist_per_length_per_cycle[int(row['len'])][int(row['cycles'])] = hist_values
        return hist_per_length_per_cycle


def get_hists_pairs_by_len(hist_per_length_per_cycle):
    hist_pairs = []
    for ms_len in [5, 10, 15, 20, 25, 30]:
        syn_hist_26 = hist_per_length_per_cycle[ms_len][26]
        syn_hist_47 = hist_per_length_per_cycle[ms_len][47]
        hist_pairs.append((ms_len, syn_hist_26, syn_hist_47))
    return hist_pairs


def get_hists_pairs(filename='hist_by_ms_len_as_0_sum.tab'):
    return get_hists_pairs_by_len(read_synthetic_data(filename))