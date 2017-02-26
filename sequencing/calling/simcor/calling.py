from sequencing.calling.simcor.models_common import BestCorrelationCalledAlleles
from sequencing.calling.hist import Histogram
from sequencing.analysis.models import HistogramEntryReads
import sys
from sequencing.calling.models import MicrosatelliteAlleleSet


def get_ms_hist(dbhist, microsatellite):
    hist_dict = {}
    for her in HistogramEntryReads.objects.filter(histogram=dbhist):
        for msg in her.microsatellite_genotypes.genotypes:
            if msg.microsatellite == microsatellite:
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
    for sim_hist in sim_space:
        distance = distance_function(real_hist, sim_hist)
        if distance < min_dist:
            min_dist = distance
            best_sim_hist = sim_hist
    return best_sim_hist, min_dist


def call_microsatellite_histogram(calling_schema, dbhist, microsatellite):
    hist = get_ms_hist(dbhist, microsatellite)
    closest_sim_hist, min_dist = get_closest(hist, calling_schema.sim_hists_space, calling_schema.distance_metric)
    mas = MicrosatelliteAlleleSet.get_for_repeats(closest_sim_hist.allele_frozenset)
    bcca, created = BestCorrelationCalledAlleles.objects.get_or_create(
        histogram=dbhist,
        microsatellite=microsatellite,
        calling_scheme=calling_schema,
        defaults=dict(
            genotypes=mas,
            confidence=min_dist,
            cycle=closest_sim_hist.simulation_cycle),
    )
    return bcca
