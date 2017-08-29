from sequencing.calling.simcor.models_common import ProportionalMicrosatelliteAlleleSet, ProportionStepModelMixin
import sys
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


def get_population_kernels(genotypes, allele_number=2, minimal_distance_between_peaks=3, case=1, filter_ones=False):
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
            min_prop=0.2)
    elif case == 2:
        return better_get_far_apart_highest_peaks(
            h,
            minimal_distance_between_peaks=minimal_distance_between_peaks,
            min_prop=0.2)
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
    for tup in pairwise_overlap(peaks):
        p1, p2 = tup
        p1_max = min(p1 + max_distance_from_peak, p1 + (p2 - p1) // 2)
        p2_min = max(p1, p2 - max_distance_from_peak, p1 + (p2 - p1) // 2)
        yield p1_max
        yield p2_min
    yield p2 + max_distance_from_peak


def get_peaks_ranges(peaks, max_distance_from_peak):
    for t in pairwise_not_overlap(get_indices(peaks, max_distance_from_peak)):
        yield range(*t)


def split_genotypes(cas, max_distance_from_peak=2, case=1, filter_ones=False):
    peaks = get_population_kernels(
        cas, allele_number=2, minimal_distance_between_peaks=3, case=case, filter_ones=filter_ones)
    if peaks is None or len(peaks) != 2:
        return
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
