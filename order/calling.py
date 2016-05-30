import sys
import numpy as np
from order.utils.parsers import uncalled_inputs
from .preprocessing import flatten_index, inflate_index
from .fitting import match_cycles
from .hist import Histogram
from itertools import combinations
from collections import defaultdict
from pickle import loads
import concurrent.futures
from itertools import repeat


def load_or_create_calling(callingfile):
    try:
        print('loading existing calling')
        f = open(callingfile, 'rb').read()
        calling = loads(f)
        print('done loading existing calling')
    except:
        print('initializing new calling')
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
    :param hist:
    :param dup_sim_hist:
    :param shift_margines:
    :param max_alleles:
    :param max_distance_from_median:
    :param kwargs:
    :return:
    """
    h = Histogram(hist, **kwargs)  # normalize, nsamples, truncate, cut_peak, trim_extremes
    med = int(np.median(h.sample))
    best_score = sys.maxsize
    res = {}
    for allele_number in range(1, max_alleles+1):
        for measured_hist_shift in range(max(0, med-shift_margins), med+shift_margins):
            normalized_shifted_reads_hist = h-measured_hist_shift
            possible_hist_seeds = combinations(
                list(range(
                    max(5, measured_hist_shift-max_distance_from_median),
                    min(measured_hist_shift+max_distance_from_median, max_ms_length))), allele_number
            )
            for seeds in possible_hist_seeds:
                if measured_hist_shift == int(np.mean(seeds)):
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


def helper(tup):
    input_tup, extras = tup
    loc, cell, row_hist = input_tup
    flat_sim_hists, kwargs = extras
    sim_hists = inflate_index(flat_sim_hists)
    # print 'working on loc: {} , cell: ...{}'.format(loc, cell[-15:])
    res = call_multi_hist(row_hist, sim_hists, **kwargs)
    return loc, cell, row_hist, res


def generate_hist_calls(input_file,
                        sim_hists,
                        calling,
                        reads_threshold=50,
                        workers=1,
                        **kwargs):
    """
        max_alleles=2,
        max_distance_from_median=30,

        shift_margins=3
        nsamples=None
        method='cor',
        score_threshold=0.006,
        min_cycles=0,
        max_cycles=80,
        max_ms_length=60
        normalize=True,
        truncate=False,
        cut_peak=False,
        trim_extremes=False):
    :param input_file:
    :param sim_hists:
    :param calling:
    :param kwargs:
    :return:
    """
    flat_sim_hists = flatten_index(sim_hists)
    print('starting {} auxiliary process'.format(workers))
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for result in executor.map(helper,
                                   zip(
                                       uncalled_inputs(input_file,
                                                       calling,
                                                       reads_threshold=reads_threshold),
                                       repeat((flat_sim_hists, kwargs))
                                   )):
            loc, cell, row_hist, res = result
            yield loc, cell, row_hist, res
