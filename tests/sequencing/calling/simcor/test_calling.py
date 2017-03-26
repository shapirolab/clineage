import pytest
from targeted_enrichment.planning.models import Microsatellite
from sequencing.calling.simcor.models_common import BestCorrelationProportionalCalledAlleles, \
    ProportionalMicrosatelliteAlleleSet, BestCorrelationProportionalHighestPeakCalledAlleles


@pytest.mark.django_db
def test_mono_calling(histograms_and_calling_solutions_d, simcormonoschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = simcormonoschema.call_ms_hist(dbhist, ms)
                alleles = bcca.genotypes.alleles
                solution_alleles = set(a for a, p in proportional_solution_alleles)
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
                solution_alleles = set(a for a, p in proportional_solution_alleles)
                assert alleles == solution_alleles  # assert called alleles against solution alleles


@pytest.mark.django_db
def test_bi_proportional_calling(histograms_and_calling_solutions_d, simcorbipropschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = simcorbipropschema.call_ms_hist(dbhist, ms)
                assert isinstance(bcca, BestCorrelationProportionalCalledAlleles)
                bcca = BestCorrelationProportionalCalledAlleles.objects.get(pk=bcca.pk)  # get rid of possible overriden field
                pmas = bcca.genotypes.proportionalmicrosatellitealleleset
                assert isinstance(pmas, ProportionalMicrosatelliteAlleleSet)
                proportional_alleles = pmas.alleles
                result_allels = set(a for a, p in proportional_alleles)
                solution_alleles = set(a for a, p in proportional_solution_alleles)
                assert result_allels == solution_alleles  # assert called alleles against solution alleles
                assert proportional_alleles == proportional_solution_alleles


@pytest.mark.django_db
def test_bi_proportional_highest_peak_calling(histograms_and_calling_solutions_d, simcorbiprophighpeakschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = simcorbiprophighpeakschema.call_ms_hist(dbhist, ms)
                assert isinstance(bcca, BestCorrelationProportionalHighestPeakCalledAlleles)
                bcca = BestCorrelationProportionalHighestPeakCalledAlleles.objects.get(pk=bcca.pk)  # get rid of possible overriden field
                pmas = bcca.genotypes.proportionalmicrosatellitealleleset
                assert isinstance(pmas, ProportionalMicrosatelliteAlleleSet)
                proportional_alleles = pmas.alleles
                result_allels = set(a for a, p in proportional_alleles)
                solution_alleles = set(a for a, p in proportional_solution_alleles)
                assert result_allels == solution_alleles  # assert called alleles against solution alleles
                assert proportional_alleles == proportional_solution_alleles
