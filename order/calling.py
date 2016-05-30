import sys
import numpy as np
from order.utils.parsers import uncalled_inputs
from .preprocessing import flatten_index, inflate_index
from .fitting import match_cycles
from .hist import Histogram
from itertools import combinations
# import concurrent.futures
from itertools import repeat
from frogress import bar


def get_proportions(diff, median, length_sensitivity=0.21, diff_sensetivity=0.65, steps=100):
    assert length_sensitivity > 0
    steps_from_equal_proportions = min(int(((diff**diff_sensetivity/(length_sensitivity*median))**3)*steps), steps/2)
    proportions = []
    for i in range(steps_from_equal_proportions):
        step = i
        proportions.append((steps/2-step)/float(steps))
        proportions.append((steps/2+step)/float(steps))
    if diff == 0:
        return [0.5]
    return proportions


def float_intersect(floats_a, floats_b, epsilon=0.000001):
    """
    forgiving floats comparison. returns the floats_a version of the values
    :param floats_a:
    :param floats_b:
    :param epsilon:
    :return:
    """
    intersection = set()
    for float_a in floats_a:
        for float_b in floats_b:
            if abs(float_a - float_b) < epsilon:
                intersection.add(float_a)
    return intersection


def filter_seeds_and_proportions(dup_sim_hist_for_seed, length_sensitivity=5.0, diff_sensetivity=1.0, steps=100):
    for seeds_and_proportions in dup_sim_hist_for_seed:
        seed_a, seed_b = seeds_and_proportions
        allele_a, proportion_a = seed_a
        allele_b, proportion_b = seed_b
        diff = abs(allele_a - allele_b)
        max_len = max([allele_a, allele_b])
        allowed_proportions = get_proportions(diff, max_len, length_sensitivity=length_sensitivity, diff_sensetivity=diff_sensetivity, steps=steps)
        if float_intersect([proportion_a, proportion_b], allowed_proportions):
            yield seeds_and_proportions


def call_multi_hist(hist,
                    dup_sim_hist,
                    shift_margins=15,
                    max_alleles=2,
                    max_distance_from_median=30,
                    max_ms_length=60,
                    proportional=False,
                    length_sensitivity=5.0,
                    diff_sensetivity=1.0,
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
            for seeds in set([frozenset(possible_hist_seed) for possible_hist_seed in possible_hist_seeds]).intersection(set(dup_sim_hist.keys())):
                if measured_hist_shift == int(np.mean(list(seeds))):
                    if not proportional:
                        c, s, best_sim_hist = match_cycles(normalized_shifted_reads_hist,
                                                           dup_sim_hist[frozenset(seeds)],
                                                           reads=h.nsamples,
                                                           proportional=proportional,
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
                    else:
                        for seeds_and_proportions in filter_seeds_and_proportions(dup_sim_hist[frozenset(seeds)], length_sensitivity=length_sensitivity, diff_sensetivity=diff_sensetivity):
                            c, s, best_sim_hist = match_cycles(normalized_shifted_reads_hist,
                                                           dup_sim_hist[frozenset(seeds)][seeds_and_proportions],
                                                           reads=h.nsamples,
                                                           **kwargs)
                            if best_score > s:
                                best_score = s
                                res = {
                                    'seeds_and_proportions': seeds_and_proportions,
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
    res = call_multi_hist(row_hist, sim_hists, proportional=False, **kwargs)
    return loc, cell, row_hist, res


# def generate_hist_calls(input_file,
#                         sim_hists,
#                         calling,
#                         reads_threshold=50,
#                         workers=1,
#                         **kwargs):
#     """
#         max_alleles=2,
#         max_distance_from_median=30,
#
#         shift_margins=3
#         nsamples=None
#         method='cor',
#         score_threshold=0.006,
#         min_cycles=0,
#         max_cycles=80,
#         max_ms_length=60
#         normalize=True,
#         truncate=False,
#         cut_peak=False,
#         trim_extremes=False):
#     :param input_file:
#     :param sim_hists:
#     :param calling:
#     :param kwargs:
#     :return:
#     """
#     flat_sim_hists = flatten_index(sim_hists)
#     print 'starting {} auxiliary process'.format(workers)
#     with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
#         for result in executor.map(helper,
#                                    izip(
#                                        uncalled_inputs(input_file,
#                                                        calling,
#                                                        reads_threshold=reads_threshold),
#                                        repeat((flat_sim_hists, kwargs))
#                                    )):
#             loc, cell, row_hist, res = result
#             yield loc, cell, row_hist, res


def generate_calling_file(input_file,
                          sim_hists,
                          calling,
                          **kwargs):
    """
        workers=1,
        max_alleles=2,
        max_distance_from_median=30,
        reads_threshold=50
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
    for loc, cell, row_hist, res in bar(generate_hist_calls(input_file, sim_hists, calling, **kwargs)):
        calling[loc][cell] = res
    return calling