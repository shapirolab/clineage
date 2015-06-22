from collections import defaultdict, Counter
from order.binomial_sim import dyn_prob
from order.hist import Histogram
from itertools import combinations, product
import numpy as np
from order.optimize_probs import dyn_mat_model
from scipy.stats import binom
from frogress import bar

def generate_bin_hist_pure_optimized(d,
                           cycles,
                           ups,
                           dws,
                           margines = 20,
                           sample_depth=10000,
                           normalize=False,
                           truncate=False,
                           cut_peak=False,
                           trim_extremes=False,
                           **kwargs):
    # ups=[lambda d:0.00005*d**2 - 0.0009*d + 0.0036],
    # dws=[lambda d:0.00009*d**2 - 0.00003*d - 0.0013],
    up = ups[0]
    dw = dws[0]
    upb = binom(cycles, min(max(0, up(d)),1.0))
    dwb = binom(cycles, min(max(0, dw(d)),1.0))
    max_mean = max(upb.mean(), dwb.mean())
    try:
        bin_margines = int(round(max_mean)) + margines
    except:
        print up(d), dw(d), d, cycles, upb.mean(), dwb.mean()
        raise
    n = np.convolve(upb.pmf(range(bin_margines)), dwb.pmf(range(bin_margines))[::-1])
    nd = {i:n[i] for i in range(bin_margines*2-1)}
    nh = Histogram(nd,
                   normalize=normalize,
                   nsamples=sample_depth,
                   truncate=truncate,
                   cut_peak=cut_peak,
                   trim_extremes=trim_extremes
                   ) - (bin_margines - 1)
    nh.truncate(p=0.0001)
    nh.normalize()
    nh.clean_zero_entries()
    return nh


def generate_dyn_hist(d,
                      cycles,
                      up,
                      dw,
                      sample_depth=10000,
                      normalize=True,
                      truncate=False,
                      cut_peak=False,
                      trim_extremes=False,
                      **kwargs):
    dyn_hist = dyn_prob(cycles,
                        d,
                        up,
                        dw,
                        nsamples=sample_depth,
                        normalize=normalize,
                        truncate=truncate,
                        cut_peak=cut_peak,
                        trim_extremes=trim_extremes)
    dyn_hist.truncate(p=0.0001)
    dyn_hist.normalize()
    dyn_hist.clean_zero_entries()
    return dyn_hist - d


def generate_mat_hist(d,
                      cycles,
                      ups,
                      dws,
                      sample_depth=10000,
                      normalize=True,
                      truncate=False,
                      cut_peak=False,
                      trim_extremes=False,
                      **kwargs):
    values = dyn_mat_model(ups, dws, d, cycles)
    h = Histogram({i: values[i] for i in range(100)},
                  nsamples=sample_depth,
                  normalize=normalize,
                  truncate=truncate,
                  cut_peak=cut_peak,
                  trim_extremes=trim_extremes)
    h.truncate(p=0.0001)
    h.normalize()
    h.clean_zero_entries()
    return h - d


def get_method(method):
    if method == 'bon':
        return generate_bin_hist_pure_optimized
    if method == 'dyn':
        return generate_dyn_hist
    if method == 'mat':
        return generate_mat_hist
    print 'unknown method'
    raise


def generate_hist(d, cycles, method, **kwargs):
    generate_method_hist = get_method(method)
    return generate_method_hist(d, cycles, **kwargs)


def generate_sim_hists(max_ms_length=60,
                       max_cycles=90,
                       method='bin',
                       **kwargs):
    """
    Creates a dictionary of simulated histogram indexed by [original_length][amplification_cycle]
    :param max_ms_length:
    :param max_cycles:
    :param method:
    :param kwargs:
    :return:
    """
    sim_hists = defaultdict(dict)
    for d in bar(range(max_ms_length)):
        for cycles in range(max_cycles):
            sim_hists[d][cycles] = generate_hist(d, cycles, method, **kwargs)
    return sim_hists


def generate_duplicate_sim_hist(sim_hists, max_alleles=2):
    dup_sim_hist = defaultdict(lambda: defaultdict(dict))
    for allele_number in range(1, max_alleles+1):
        for seeds in combinations(sim_hists.keys(), allele_number):  # iterate over all possible original lengths
            shift = int(np.mean(seeds))
            for cycles in bar(sim_hists[0].keys()):
                first_seed = seeds[0]  # initial microsatellite length (length = seed[0]
                sum_hist = sim_hists[first_seed][cycles] + first_seed  # add generated histogram based on seed and shift it on the x axis according to the seed's value (simulated microsatellite length
                for seed in seeds[1:]:
                    sum_hist = sum_hist + (sim_hists[seed][cycles] + seed)  # iteratively add generated histograms shifted acoording to their simulated seed
                dup_sim_hist[frozenset(seeds)][cycles] = sum_hist - shift  # shift the sum of all histograms so that thier based around zero. The shift
    return dup_sim_hist


def generate_simulated_proportional_alleles(seeds, cycles, proprtions, method, **kwargs):
    assert sum(proprtions) == 1
    assert len(seeds) == len(cycles) == len(proprtions)
    signals = []
    for seed, cycles, proportion in zip(seeds, cycles, proprtions):
        hist = generate_hist(seed, cycles, method, **kwargs).ymul(proportion) + seed
        signals.append(hist)
    hist_sum = signals[0]
    for signal in signals[1:]:
        hist_sum = hist_sum.asym_add(signal)
    hist_sum.normalize()
    shift = int(np.mean(seeds))
    return hist_sum - shift


def generate_biallelic_reads_of_multiple_proportions(max_ms_length=60, max_cycles=90, method='bin', **kwargs):
    for seeds in combinations(range(max_ms_length), 2):
        for cycles in range(max_cycles):
            for p1 in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
                dup_sim_hist[frozenset(seeds)][cycles][p1] = generate_simulated_proportional_alleles(seeds, cycles, [p1, 1-p1], method, **kwargs)

def generate_sim_hists_of_up_to_k_alleles(**kwargs):
    """
        method='bin'
        max_ms_length=60,
        max_cycles=90,
        ups=[lambda x: 0.003],
        dws=[lambda x: 0.022],
        sample_depth=10000,
        normalize=True,
        truncate=False,
        cut_peak=False,
        trim_extremes=False
        max_alleles = 2
    :param kwargs:
    :return:
    """
    sim_hists = generate_sim_hists(**kwargs)
    dup_sim_hist = generate_duplicate_sim_hist(sim_hists, kwargs['max_alleles'])
    return dup_sim_hist


def flatten_index(sim_hists):
    flat_dict = {}
    for root in sim_hists:
        for cycle in sim_hists[root]:
            flat_dict[(root, cycle)] = sim_hists[root][cycle]
    return flat_dict


def inflate_index(sim_hists):
    inflated_dict = defaultdict(lambda: defaultdict(dict))
    for root, cycle in sim_hists:
        inflated_dict[root][cycle] = sim_hists[(root, cycle)]
    return inflated_dict
