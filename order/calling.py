import sys
import numpy as np

from fitting import match_cycles
from hist import Histogram


def call_hist(hist, sim_hists, min_cycles=0, max_cycles=50, method='cor', shift_margines=3, normalize=True, nsamples=None, trunc=False, cut_peak=False, trim_extremes=False):
    h = Histogram(hist, normalize=normalize, nsamples=nsamples, trunc=trunc, cut_peak=cut_peak, trim_extremes=trim_extremes)
    if h.nsamples == 0:
        print loc, cell
        raise
    med = int(np.median(h.sample))
    best_score = sys.maxint
    res = {}
    for d in range(max(0, med-shift_margines), med+shift_margines):
        normalized_shifted_reads_hist = h-d
        c, s, best_sim_hist = match_cycles(normalized_shifted_reads_hist, sim_hists, d, method=method, reads=h.nsamples, min_cycles=min_cycles, max_cycles=max_cycles)
        #print med, d, c, s
        if best_score > s:
            best_score = s
            res = {
                  'shift':d,
                  'cycle':c,
                  'score':s,
                  'median':med,
                  'reads':h.nsamples
                  }
    return res