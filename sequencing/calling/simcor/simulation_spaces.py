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


def remove_points_close_to_top(hs, number_of_points, distance=1):
    """
    Removes points on the histogram if they are withing the given range, and returns the histogram without them
    Searches if there's a point in the given range, and removes it,
    :param hs: sorted histogram
    :param keep_point_position: the position point you wish to keep
    :param distance: the distance to remove
    :return: remaining list
    """
    for steps_from_point in range(1, distance):
        if [item for item in hs if item[0] == hs[number_of_points][0] + steps_from_point]:
            hs.remove([item for item in hs if item[0] == hs[number_of_points][0] + steps_from_point][0])
        if [item for item in hs if item[0] == hs[number_of_points][0] - steps_from_point]:
            hs.remove([item for item in hs if item[0] == hs[number_of_points][0] - steps_from_point][0])
    return hs


def get_far_apart_highest_peaks(hist, allele_number=1, minimal_distance_between_peaks=1):
    """
    Identify the k highest peaks that satisfy minimal distance
    Args:
        hist: histogram
        allele_number: allele_number
        minimal_distance_between_peaks: minimum distance between the allele
    """
    hs = sorted(hist._hist.items(), key=lambda hkey: hkey[1], reverse=True)
    for allele in range(allele_number):
        if allele >= len(hs):
            break
        hs = remove_points_close_to_top(hs, allele, minimal_distance_between_peaks)
    seeds = [x for x, y in hs[:allele_number]]
    return seeds


def seeds_search_range(peaks, max_distance_between_peaks, max_ms_length ):
    search_range = itertools.product(
            *[
                range(
                    max(1, peak-max_distance_between_peaks),
                    min(max_ms_length, peak+max_distance_between_peaks+1)
                ) for peak in peaks
            ])
    yield from search_range
