import pytest
from targeted_enrichment.planning.models import Microsatellite


@pytest.mark.django_db
def test_mono_calling(histograms_and_calling_solutions_d, simcormonoschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = simcormonoschema.call_ms_hist(dbhist, ms)
                alleles = bcca.genotypes.alleles
                solution_alleles = set(a for p, a in proportional_solution_alleles)
                assert len(alleles) == 1  # Assert mono solution
                if len(solution_alleles) == 1:  # Monoallele test
                    assert alleles == solution_alleles  # assert called alleles against solution alleles
                else:  # Multiallele test
                    assert alleles & solution_alleles  # assert called alleles against solution alleles


@pytest.mark.django_db
def test_bi_calling(histograms_and_calling_solutions_d, simcorbinschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = simcorbinschema.call_ms_hist(dbhist, ms)
                alleles = bcca.genotypes.alleles
                solution_alleles = set(a for p, a in proportional_solution_alleles)
                assert alleles == solution_alleles  # assert called alleles against solution alleles
