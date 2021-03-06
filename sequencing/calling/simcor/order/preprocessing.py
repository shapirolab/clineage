from collections import defaultdict, Counter
from sequencing.calling.hist import Histogram
from itertools import combinations, combinations_with_replacement, product
import numpy as np
from frogress import bar


def generate_mat_hist(ms_len, 
        mat,
        sample_depth=10000,
        normalize=True,
        truncate=False,
        cut_peak=False,
        trim_extremes=False,
        n=50,
        **kwargs):
    h = Histogram({i: mat[ms_len, i] for i in range(n)},
                  nsamples=sample_depth,
                  normalize=normalize,
                  truncate=truncate,
                  cut_peak=cut_peak,
                  trim_extremes=trim_extremes)
    h.truncate(p=0.0001)
    h.normalize()
    h.clean_zero_entries()
    return h - ms_len


# def generate_sim_hists(max_ms_length=60,
#                        max_cycles=90,
#                        method='mat',
#                        **kwargs):
#     """
#     Creates a dictionary of simulated histogram indexed by [original_length][amplification_cycle]
#     :param max_ms_length:
#     :param max_cycles:
#     :param method:
#     :param kwargs:
#     :return:
#     """
#     sim_hists = defaultdict(dict)
#     for d in bar(range(max_ms_length)):
#         for cycles in range(max_cycles):
#             sim_hists[d][cycles] = generate_hist(d, cycles, method, **kwargs)
#     return sim_hists
#
#
# def generate_duplicate_sim_hist(sim_hists, max_alleles=2):
#     dup_sim_hist = defaultdict(lambda: defaultdict(dict))
#     for allele_number in range(1, max_alleles+1):
#         for seeds in combinations(sim_hists.keys(), allele_number):  # iterate over all possible original lengths
#             shift = int(np.mean(seeds))
#             for cycles in bar(sim_hists[0].keys()):
#                 first_seed = seeds[0]  # initial microsatellite length (length = seed[0]
#                 sum_hist = sim_hists[first_seed][cycles] + first_seed  # add generated histogram based on seed and shift it on the x axis according to the seed's value (simulated microsatellite length
#                 for seed in seeds[1:]:
#                     sum_hist = sum_hist + (sim_hists[seed][cycles] + seed)  # iteratively add generated histograms shifted acoording to their simulated seed
#                 dup_sim_hist[frozenset(seeds)][cycles] = sum_hist - shift  # shift the sum of all histograms so that thier based around zero. The shift
#     return dup_sim_hist
#
#
# def generate_simulated_proportional_alleles(seeds, cycles, proprtions, method, **kwargs):
#     assert sum(proprtions) == 1
#     assert len(seeds) == len(cycles) == len(proprtions)
#     signals = []
#     for seed, cycles, proportion in zip(seeds, cycles, proprtions):
#         hist = generate_hist(seed, cycles, method, **kwargs).ymul(proportion) + seed
#         signals.append(hist)
#     hist_sum = signals[0]
#     for signal in signals[1:]:
#         hist_sum = hist_sum.asym_add(signal)
#     hist_sum.normalize()
#     shift = int(np.mean(seeds))
#     return hist_sum - shift
#
#
# def generate_simulated_proportional_alleles_precalculated(seeds, seeds_hists, cycles, proprtions):
#     assert sum(proprtions) == 1
#     assert len(seeds) == len(cycles) == len(proprtions)
#     signals = []
#     for seed, cycles, proportion in zip(seeds, cycles, proprtions):
#         if proportion == 0.0:
#             continue
#         hist = seeds_hists[seed][cycles].ymul(proportion) + seed
#         signals.append(hist)
#     hist_sum = signals[0]
#     for signal in signals[1:]:
#         hist_sum = hist_sum.asym_add(signal)
#     hist_sum.normalize()
#     shift = int(np.mean(list(seeds)))
#     return hist_sum - shift
#
#
# def generate_biallelic_reads_of_multiple_proportions(min_cycles=20, max_cycles=90, min_ms_length=5, max_ms_length=60, method='bin', steps=100,  **kwargs):
#     # print min_ms_length, max_ms_length, min_cycles, max_cycles
#     sim_hists = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
#     for seeds in bar(list(combinations_with_replacement(list(range(min_ms_length, max_ms_length)), 2))):
#         for cycles in range(min_cycles, max_cycles):
#             for p1 in [float(x)/steps for x in range(1, steps-1)]:
#                 # print p1
#                 sim_hists[frozenset(seeds)][tuple(zip(seeds, (p1, 1 - p1)))][cycles] = generate_simulated_proportional_alleles(seeds, (cycles, cycles), (p1, 1 - p1), method, **kwargs)
#     return sim_hists
#
#
# def generate_sim_hists_of_up_to_k_alleles(**kwargs):
#     """
#         method='bin'
#         max_ms_length=60,
#         max_cycles=90,
#         ups=[lambda x: 0.003],
#         dws=[lambda x: 0.022],
#         sample_depth=10000,
#         normalize=True,
#         truncate=False,
#         cut_peak=False,
#         trim_extremes=False
#         max_alleles = 2
#     :param kwargs:
#     :return:
#     """
#     sim_hists = generate_sim_hists(**kwargs)
#     dup_sim_hist = generate_duplicate_sim_hist(sim_hists, kwargs['max_alleles'])
#     return dup_sim_hist
#
#
# def flatten_index(sim_hists):
#     flat_dict = {}
#     for root in sim_hists:
#         for cycle in sim_hists[root]:
#             flat_dict[(root, cycle)] = sim_hists[root][cycle]
#     return flat_dict
#
#
# def inflate_index(sim_hists):
#     inflated_dict = defaultdict(lambda: defaultdict(dict))
#     for root, cycle in sim_hists:
#         inflated_dict[root][cycle] = sim_hists[(root, cycle)]
#     return inflated_dict
