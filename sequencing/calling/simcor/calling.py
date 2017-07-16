from sequencing.calling.simcor.models_common import ProportionalMicrosatelliteAlleleSet, ProportionStepModelMixin, \
    Proportions
from sequencing.calling.hist import Histogram
from sequencing.analysis.models import HistogramEntryReads
import sys
from sequencing.calling.models import MicrosatelliteAlleleSet
from misc.utils import get_get_or_create


def get_ms_hist(dbhist, microsatellite):
    hist_dict = {}
    for her in HistogramEntryReads.objects.filter(histogram=dbhist):
        for msg in her.microsatellite_genotypes.genotypes:
            if msg.microsatellite == microsatellite:
                if msg.repeat_number in hist_dict:
                    hist_dict[msg.repeat_number] += her.num_reads
                else:
                    hist_dict[msg.repeat_number] = her.num_reads
    return Histogram(hist_dict)


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
    return get_get_or_create(inner, calling_schema.called_allele_class,
         histogram=dbhist,
         microsatellite=microsatellite,
         calling_scheme=calling_schema,
    )
