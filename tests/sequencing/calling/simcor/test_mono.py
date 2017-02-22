import pytest
from tests.sequencing.calling.conftest import *
from sequencing.calling.simcor.calling import call_microsatellite_histogram
from tests.sequencing.calling.simcor.conftest import *


@pytest.mark.django_db
def test_genotype_mapping(histograms_and_calling_solutions, simcor):
    for solution_alleles, dbhist, ms in histograms_and_calling_solutions:
        # do calling (histogram, microsatellite)
        # get CalledAlleles
        best_sim_dict, min_dist = call_microsatellite_histogram(dbhist, ms)
        allele = best_sim_dict
        assert ca.genotypes.alleles == set(a for p, a in solution_alleles)  # assert called alleles against solution alleles


def call_microsatellite_histogram(dbhist, microsatellite, sim_hists_space, simcor):
    hist = get_ms_hist(dbhist, microsatellite)
    best_sim_dict, min_dist = get_closest(hist, sim_hists_space(sim_hists_space), pop_dist_corr_numpy)
    mas = MicrosatelliteAlleleSet.get_for_repeats(best_sim_dict['ms_len'])
    ca, created = CalledAlleles.objects.get_or_create(
        histogram=dbhist,
        microsatellite=microsatellite,
        calling_scheme=simcor,
        defaults=dict(genotype=mas),
    )
    return best_sim_dict, min_dist
