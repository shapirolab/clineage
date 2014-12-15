import sys
import numpy as np

from fitting import match_cycles
from hist import Histogram
from itertools import combinations
from collections import defaultdict
from cloud.serialization.cloudpickle import loads


def load_or_create_calling(callingfile):
    try:
        f = open(callingfile,'rb').read()
        calling = loads(f)
    except:
        calling = defaultdict(lambda: defaultdict(dict))
    return calling


def call_multi_hist(hist,
                  dup_sim_hist,
                  shift_margins=15,
                  max_alleles=2,
                  max_distance_from_median=30,
                  max_ms_length=60,
                  **kwargs
                  ):
    """
        normalize=True,
        nsamples=None,
        truncate=False,
        cut_peak=False,
        trim_extremes=False,
        min_cycles=0,
        max_cycles=50,
        method='cor',
        max_ms_length=60
    :param hist:
    :param dup_sim_hist:
    :param shift_margines:
    :param max_alleles:
    :param max_distance_from_median:
    :param kwargs:
    :return:
    """
    h = Histogram(hist, **kwargs) #normalize, nsamples, truncate, cut_peak, trim_extremes
    med = int(np.median(h.sample))
    best_score = sys.maxint
    res = {}
    for allele_number in range(1, max_alleles+1):
        for measured_hist_shift in range(max(0, med-shift_margins), med+shift_margins):
            normalized_shifted_reads_hist = h-measured_hist_shift
            possible_hist_seeds = combinations(
                range(
                    max(0, measured_hist_shift-max_distance_from_median),
                    min(measured_hist_shift+max_distance_from_median, max_ms_length)
                ), allele_number
            )
            for seeds in possible_hist_seeds:
                c, s, best_sim_hist = match_cycles(normalized_shifted_reads_hist,
                                                   dup_sim_hist[frozenset(seeds)],
                                                   reads=h.nsamples,
                                                   **kwargs)
                if best_score > s:
                    best_score = s
                    res = {
                          'shifts': seeds,
                          'cycle': c,
                          'score': s,
                          'median': med,
                          'reads': h.nsamples
                          }
    return res