from sequencing.calling.simcor.models_common import ProportionalMicrosatelliteAlleleSet, ProportionStepModelMixin
from sequencing.calling.hist_dist import derived_proportions_dot
import sys
import numpy as np
from sequencing.calling.models import MicrosatelliteAlleleSet
from misc.utils import get_get_or_create
from collections import Counter
from sequencing.analysis.models import HistogramEntryReads, Histogram
from targeted_enrichment.amplicons.models import Amplicon
from itertools import tee

from sequencing.calling.hist import Histogram as dHistogram
from sequencing.calling.simcor.hist_analysis import get_far_apart_highest_peaks, better_get_far_apart_highest_peaks,\
    better_get_far_apart_highest_peaks_that_doesnt_hang


def get_ms_hist(dbhist, microsatellite):
    hist_dict = {}
    for her in HistogramEntryReads.objects.filter(histogram=dbhist):
        for msg in her.microsatellite_genotypes.genotypes:
            if msg.microsatellite == microsatellite:
                if msg.repeat_number in hist_dict:
                    hist_dict[msg.repeat_number] += her.num_reads
                else:
                    hist_dict[msg.repeat_number] = her.num_reads
    return dHistogram(hist_dict)


def get_closest(real_hist, sim_space, distance_function):
    """
    Measure a histogram against a simulation space and return the closest point in space
    Args:
        real_hist: Histogram object
        sim_space: SimulatedHistograms generator
        distance_function: lower-is-closer distance function

    Returns:
        SimulatedHistogram with minimal distance to real_hist
    """
    min_dist = sys.maxsize
    best_sim_hist = None
    examined = set()
    for sim_hist in sim_space:
        if sim_hist.identity in examined:
            continue
        distance = distance_function(real_hist, sim_hist)
        if distance < min_dist:
            min_dist = distance
            best_sim_hist = sim_hist
        examined.add(sim_hist.identity)
    return best_sim_hist, min_dist


def real_hist_to_rht(real_hist):
    rhv = real_hist._vec
    rhmr = np.matrix(rhv)
    rhmc = np.matrix(rhv).transpose()
    rhv_n = np.linalg.norm(rhv, ord=2)
    rht = (rhmr, rhmc, rhv_n)
    return rht

from sequencing.calling.simcor.range import contains_excluded_proportions_wrapper
def get_closest_vec_opt(real_hist, sim_space, distance_function, length_sensitivity=0.21, diff_sensetivity=0.65):
    """
    Measure a histogram against a simulation space and return the closest point in space
    Args:
        real_hist: Histogram object
        sim_space: SimulatedHistograms generator
        distance_function: lower-is-closer distance function

    Returns:
        SimulatedHistogram with minimal distance to real_hist
    """
    assert distance_function==derived_proportions_dot
    min_dist = sys.maxsize
    best_sim_hist = None
    examined = set()
    rht = real_hist_to_rht(real_hist)
    for sim_hist in sim_space:
        if sim_hist.identity in examined:
            continue
        conf, p = distance_function(rht, sim_hist)
        if conf is None:
            continue
        distance = 1 - conf
        if distance < min_dist:
            # (frozenset({(6, Decimal('1.00000')), (5, Decimal('0.00000'))}), cyc)
            if contains_excluded_proportions_wrapper(sim_hist.identity, length_sensitivity=length_sensitivity, diff_sensetivity=diff_sensetivity):
                continue
            min_dist = distance
            best_sim_hist = sim_hist
            best_sim_hist._p1 = p
            best_sim_hist._p2 = 1 - p
            best_sim_hist._ms_lens_and_proportions = (best_sim_hist._a1, best_sim_hist._p1), (best_sim_hist._a2, best_sim_hist._p2)
            best_sim_hist._alleles_to_proportions = {a: p for a, p in best_sim_hist._ms_lens_and_proportions if p > 0}
        examined.add(sim_hist.identity)
    return best_sim_hist, min_dist


def get_closest_vec_opt_mms(real_hist, sim_space, distance_function, length_sensitivity=0.21, diff_sensetivity=0.65, eps=0.05):
    """
    Measure a histogram against a simulation space and return the closest point in space
    Args:
        real_hist: Histogram object
        sim_space: SimulatedHistograms generator
        distance_function: lower-is-closer distance function

    Returns:
        SimulatedHistogram with minimal distance to real_hist
    """
    assert distance_function==derived_proportions_dot
    min_dist = sys.maxsize
    best_sim_hist = None
    second_min_dist = sys.maxsize
    second_best_sim_hist = None
    examined = set()
    rht = real_hist_to_rht(real_hist)
    for sim_hist in sim_space:
        if sim_hist.identity in examined:
            continue
        conf, p = distance_function(rht, sim_hist, eps=eps)
        if conf is None:
            continue
        sim_hist._p1 = p
        sim_hist._p2 = 1 - p
        sim_hist._ms_lens_and_proportions = (sim_hist._a1, sim_hist._p1), (
            sim_hist._a2, sim_hist._p2)
        sim_hist._alleles_to_proportions = {a: p for a, p in sim_hist._ms_lens_and_proportions if p > 0}

        distance = 1 - conf
        if distance < min_dist:
            # (frozenset({(6, Decimal('1.00000')), (5, Decimal('0.00000'))}), cyc)
            if contains_excluded_proportions_wrapper(sim_hist.identity, length_sensitivity=length_sensitivity,
                                                     diff_sensetivity=diff_sensetivity):
                continue
            if best_sim_hist is not None:
                if sim_hist.allele_frozenset != best_sim_hist.allele_frozenset:
                    second_min_dist = min_dist
                    second_best_sim_hist = best_sim_hist
            min_dist = distance
            best_sim_hist = sim_hist
        elif second_min_dist is None or distance < second_min_dist:
            if contains_excluded_proportions_wrapper(sim_hist.identity, length_sensitivity=length_sensitivity,
                                                     diff_sensetivity=diff_sensetivity):
                continue
            if sim_hist.allele_frozenset != best_sim_hist.allele_frozenset:
                second_min_dist = distance
                second_best_sim_hist = sim_hist
        examined.add(sim_hist.identity)
    return (best_sim_hist, min_dist), (second_best_sim_hist, second_min_dist)


def call_microsatellite_histogram(calling_schema, dbhist, microsatellite):
    def inner(raise_or_create_with_defaults):
        hist = get_ms_hist(dbhist, microsatellite)
        closest_sim_hist, min_dist = calling_schema.find_best_in_space(hist)
        mas = MicrosatelliteAlleleSet.get_for_alleles(closest_sim_hist.allele_frozenset)
        if isinstance(calling_schema, ProportionStepModelMixin):
            mas = ProportionalMicrosatelliteAlleleSet.get_for_proportional_alleles(mas, closest_sim_hist.alleles_to_proportions)
        return raise_or_create_with_defaults(
            genotypes=mas,
            confidence=min_dist,
            cycle=closest_sim_hist.simulation_cycle,
        )
    return get_get_or_create(inner,
                             calling_schema.called_allele_class,
                             histogram=dbhist,
                             microsatellite=microsatellite,
                             calling_scheme=calling_schema,
                             )


def get_amplicon_by_ms(ms):
    return Amplicon.objects.filter(slice__contains=ms.slice)


def get_amplicons_by_sr(sr):
    return set(amp.id for amp in sr.library.subclass.amplicons)


def get_ms_amplicon(ms, sr_amps):
    amplicons = get_amplicon_by_ms(ms)
    actual_amp = set(amp.id for amp in amplicons) & sr_amps
    ampid = actual_amp.pop()
    return Amplicon.objects.select_subclasses().get(id=ampid)


def get_population_kernels(genotypes, allele_number=2, minimal_distance_between_peaks=3, case=1, filter_ones=False,
                           min_prop=0.2):
    if filter_ones:
        # Ignore single occurances
        h = dHistogram(
            {k: v for k, v in Counter([a for ca in genotypes for a in ca.genotypes.alleles]).items() if v > 1}
        )
        if len(h.keys()) == 0:
            return None
    else:
        h = dHistogram(Counter([a for ca in genotypes for a in ca.genotypes.alleles]))

    if case == 1:
        return get_far_apart_highest_peaks(
            h,
            allele_number=allele_number,
            minimal_distance_between_peaks=minimal_distance_between_peaks,
            min_prop=min_prop)
    elif case == 2:
        return better_get_far_apart_highest_peaks(
            h,
            minimal_distance_between_peaks=minimal_distance_between_peaks,
            min_prop=min_prop)
    elif case == 3:
        return better_get_far_apart_highest_peaks_that_doesnt_hang(
            h,
            minimal_distance_between_peaks=minimal_distance_between_peaks)


def pairwise_overlap(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def pairwise_not_overlap(iterable):
    """s -> (s0, s1), (s2, s3), (s4, s5), ..."""
    a = iter(iterable)
    return zip(a, a)


def get_indices(peaks, max_distance_from_peak):
    """This method gets the peaks ands returns ranges of proper values near them"""
    p1 = peaks[0]
    yield max(0, p1 - max_distance_from_peak)
    if len(peaks) > 1:
        for tup in pairwise_overlap(peaks):
            p1, p2 = tup
            p1_max = min(p1 + max_distance_from_peak, p1 + (p2 - p1) // 2)
            p2_min = max(p1, p2 - max_distance_from_peak, p1 + (p2 - p1) // 2)
            yield p1_max
            yield p2_min
        yield p2 + max_distance_from_peak
    else:
        yield p1 + max_distance_from_peak


def get_peaks_ranges(peaks, max_distance_from_peak):
    for t in pairwise_not_overlap(get_indices(peaks, max_distance_from_peak)):
        yield range(*t)


def split_genotypes(cas, max_distance_from_peak=2, case=1, filter_ones=False, min_prop=0.2, filter_single_peak=True, minimal_distance_between_peaks=3):
    peaks = get_population_kernels(
        cas, allele_number=2, minimal_distance_between_peaks=minimal_distance_between_peaks, case=case, filter_ones=filter_ones, min_prop=min_prop)
    if peaks is None or len(peaks) > 2:
        return None
    if len(peaks) == 1 and filter_single_peak:
        return None
    peaks.sort()
    peaks_by_range = {p: prange for p, prange in zip(peaks, get_peaks_ranges(peaks, max_distance_from_peak))}
    calling_assignments = dict()
    for ca in cas:
        assigned_alleles = dict()
        for a in ca.genotypes.alleles:
            for p in peaks_by_range:
                if a in peaks_by_range[p]:
                    assigned_alleles[a] = p
                    break
            else:
                pass  # TODO: allele was not assigned to a window, consider exception
        calling_assignments[ca] = assigned_alleles
    return calling_assignments
