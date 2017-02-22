from sequencing.calling.simcor.models import BestCorrelationCalledAlleles
from sequencing.calling.hist import Histogram
from sequencing.analysis.models import HistogramEntryReads
import sys
from sequencing.calling.hist_dist import pop_dist_corr_numpy
from sequencing.calling.models import MicrosatelliteAlleleSet
from sequencing.calling.simcor.simulation_spaces import mono_sim_hists_space_generator, bi_sim_hists_space_generator,\
    proportional_bi_sim_hists_space_generator


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


def call_microsatellite_histogram(dbhist, microsatellite, calling_schema):
    """
    Mono
    Args:
        dbhist:
        microsatellite:
        calling_schema:

    Returns:

    """
    hist = get_ms_hist(dbhist, microsatellite)
    sim_hists_space = mono_sim_hists_space_generator(
        calling_schema.simulations,
        calling_schema.seeds_and_cycles)
    closest_sim_hist, min_dist = get_closest(hist, sim_hists_space, pop_dist_corr_numpy)
    mas = MicrosatelliteAlleleSet.get_for_repeats([closest_sim_hist.ms_len])
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


def call_microsatellite_histogram_symetric_bi(dbhist, microsatellite, calling_schema):
    """
    Bi without proportions
    Args:
        dbhist:
        microsatellite:
        calling_schema:

    Returns:

    """
    hist = get_ms_hist(dbhist, microsatellite)
    sim_hists_space = bi_sim_hists_space_generator(
        calling_schema.simulations,
        calling_schema.seeds_and_cycles)
    closest_sim_hist, min_dist = get_closest(hist, sim_hists_space, pop_dist_corr_numpy)
    mas = MicrosatelliteAlleleSet.get_for_repeats(closest_sim_hist.ms_lens)
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