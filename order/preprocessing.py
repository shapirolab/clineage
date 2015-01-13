from collections import defaultdict, Counter
from binomial_sim import sim, dyn_prob
from frogress import bar as tqdm
from order.hist import Histogram
from itertools import combinations
import numpy as np


def generate_bin_hist(d,
                  cycles,
                  up=lambda x: 0.003,
                  dw=lambda x: 0.022,
                  sample_depth=10000,
                  normalize=True,
                  truncate=False,
                  cut_peak=False,
                  trim_extremes=False,
                  **kwargs):
    up_p = up(d)
    dw_p = dw(d)
    z = sim(cycles, up_p, dw_p)
    return Histogram(
        Counter(z.rand(sample_depth)),
        normalize=normalize,
        nsamples=sample_depth,
        truncate=truncate,
        cut_peak=cut_peak,
        trim_extremes=trim_extremes
    )


def generate_dyn_hist(d,
                  cycles,
                  up=lambda x: 0.003,
                  dw=lambda x: 0.022,
                  sample_depth=10000,
                  normalize=True,
                  truncate=False,
                  cut_peak=False,
                  trim_extremes=False,
                  **kwargs):
    dyn_hist = dyn_prob(cycles, d, up, dw)
    dyn_hist.normalize()
    return dyn_hist - d


def get_method(method):
    if method == 'bin':
        return generate_bin_hist
    if method == 'dyn':
        return generate_dyn_hist
    return


def generate_sim_hists(max_ms_length=60,
                           max_cycles=90,
                           method='bin',
                           **kwargs):
    generate_hist = get_method(method)
    sim_hists = defaultdict(dict)
    for d in tqdm(range(max_ms_length)):
        for cycles in range(max_cycles):
            sim_hists[d][cycles] = generate_hist(d, cycles, **kwargs)
    return sim_hists


def generate_sim_hists(method='bin',
                       max_ms_length=60,
                       max_cycles=90,
                       up=lambda x: 0.003,
                       dw=lambda x: 0.022,
                       sample_depth=10000,
                       normalize=True,
                       truncate=False,
                       cut_peak=False,
                       trim_extremes=False,
                       **kwargs):
    """
    """
    return generate_sim_hists(method=method,
                              max_ms_length=max_ms_length,
                              max_cycles=max_cycles,
                              up=up,
                              dw=dw,
                              sample_depth=sample_depth,
                              normalize=normalize,
                              truncate=truncate,
                              cut_peak=cut_peak,
                              trim_extremes=trim_extremes)


def generate_duplicate_sim_hist(sim_hists, max_alleles=2):
    dup_sim_hist = defaultdict(lambda: defaultdict(dict))
    for allele_number in range(1, max_alleles+1):
        for seeds in combinations(sim_hists.keys(), allele_number):
            shift = int(np.mean(seeds))
            for cycles in tqdm(sim_hists[0].keys()):#TODO: get sim_hists parameters in call?
                first_seed = seeds[0]
                sum_hist = sim_hists[first_seed][cycles] + first_seed
                for seed in seeds[1:]:
                    sum_hist = sum_hist + (sim_hists[seed][cycles] + seed)
                dup_sim_hist[frozenset(seeds)][cycles] = sum_hist - shift
    return dup_sim_hist


def generate_sim_hists_of_up_to_k_alleles(**kwargs):
    """
        method='bin'
        max_ms_length=60,
        max_cycles=90,
        up=lambda x: 0.003,
        dw=lambda x: 0.022,
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