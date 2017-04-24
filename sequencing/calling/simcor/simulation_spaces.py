import itertools

from sequencing.calling.hist import Histogram
from sequencing.calling.multi_hists import MonoSimulatedHistogram, MultiSimulatedHistogram, \
    ProportionalMultiSimulatedHistogram


def mono_sim_hists_space_generator(sim_hists_dict, seeds_and_cycles):
    """
    Generates simulated histograms with reference MS length and simulation cycles
    Args:
        sim_by_cyc: SimultaionsByCycles class instance
        seeds_and_cycles: a generator for the desired seeds and cycles that will be simulated
            [(syn_len, sim_cyc), (syn_len, sim_cyc), ...]
    """
    for syn_len, sim_cyc in seeds_and_cycles:
        yield MonoSimulatedHistogram(
            ms_len=syn_len,
            simulation_cycle=sim_cyc,
            simulated_hist=sim_hists_dict[syn_len][sim_cyc]
        )


def bi_sim_hists_space_generator(sim_hists_dict, seeds_and_cycles):
    """
    Generates simulated histograms with reference MS length and simulation cycles
    Args:
        sim_by_cyc: SimultaionsByCycles class instance
        seeds_and_cycles: a generator for the desired seeds and cycles that will be simulated
            [(frozenset({syn_len, syn_len}), sim_cyc), (frozenset({syn_len, syn_len}), sim_cyc), ...]
    """
    for syn_seeds, sim_cyc in seeds_and_cycles:
        yield MultiSimulatedHistogram(
            ms_lens=syn_seeds,
            simulation_cycle=sim_cyc,
            simulated_hist=sum(
                Histogram(
                    sim_hists_dict[syn_len][sim_cyc]
                ) for syn_len in syn_seeds)
        )


def proportional_bi_sim_hists_space_generator(sim_hists_dict, seeds_and_cycles):
    """
    Generates simulated histograms with reference MS length and simulation cycles
    Args:
        sim_by_cyc: SimultaionsByCycles class instance
        seeds_and_cycles: a generator for the desired seeds and cycles that will be simulated
            [(frozenset({(syn_len, p), (syn_len, p)}), sim_cyc), ...]
    """
    for ms_lens_and_proportions, sim_cyc in seeds_and_cycles:
        model_hist = Histogram(dict())
        for syn_len, p in ms_lens_and_proportions:
            model_hist = model_hist.asym_add(sim_hists_dict[syn_len][sim_cyc].ymul(p))
        yield ProportionalMultiSimulatedHistogram(
            ms_lens_and_proportions=ms_lens_and_proportions,
            simulation_cycle=sim_cyc,
            simulated_hist=model_hist,
        )


def seeds_search_range(peaks, max_distance_between_peaks, max_ms_length ):
    search_range = itertools.product(
            *[
                range(
                    max(1, peak-max_distance_between_peaks),
                    min(max_ms_length, peak+max_distance_between_peaks+1)
                ) for peak in peaks
            ])
    yield from search_range
