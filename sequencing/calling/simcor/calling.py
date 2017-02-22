from sequencing.calling.simcor.models import BestCor
from sequencing.calling.hist import Histogram
from sequencing.analysis.models import HistogramEntryReads
import sys
from sequencing.calling.hist_dist import pop_dist_corr_numpy
from sequencing.calling.models import MicrosatelliteAlleleSet, CalledAlleles


def get_ms_hist(dbhist, microsatellite):
    hist_dict = {}
    for her in HistogramEntryReads.objects.filter(histogram=dbhist):
        for msg in her.microsatellite_genotypes.genotypes:
            if msg.microsatellite == microsatellite:
                hist_dict[msg.repeat_number] = her.num_reads
    return Histogram(hist_dict)


def sim_hists_space_generator(sim_by_cyc):  # brute force over all mono-allelic options
    sim_hists = sim_by_cyc.get_simulations_dict()
    for syn_len in sim_hists:
        for sim_cyc in sim_hists[syn_len]:
            yield {
                'ms_len': syn_len,
                'simulation_cycle': sim_cyc,
                'simulated_hist': sim_hists[syn_len][sim_cyc]
            }


def get_closest(real_hist, sim_space, distance_function):
    min_dist = sys.maxsize
    best_sim_dict = dict()
    for sim_dict in sim_space:
        sim_hist = sim_dict['simulated_hist']
        distance = distance_function(real_hist, sim_hist)
        if distance < min_dist:
            min_dist = distance
            best_sim_dict = sim_dict
    return best_sim_dict, min_dist


# def call_microsatellite_histogram(dbhist, microsatellite, sim_hists_space):
#     hist = get_ms_hist(dbhist, microsatellite)
#     best_sim_dict, min_dist = get_closest(hist, sim_hists_space, pop_dist_corr_numpy)
#     mas = MicrosatelliteAlleleSet.get_for_repeats(best_sim_dict['ms_len'])
#     ca, created = CalledAlleles.objects.get_or_create(
#         histogram=dbhist,
#         microsatellite=microsatellite,
#         calling_scheme=,
#         defaults=dict(genotype=mas),
#     )
#     return best_sim_dict, min_dist
