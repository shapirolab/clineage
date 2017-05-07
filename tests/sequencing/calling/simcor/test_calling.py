import pytest
import decimal
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
def test_mono_calling_highest_peak(histograms_and_calling_solutions_d, simcormonoprophighpeakschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = simcormonoprophighpeakschema.call_ms_hist(dbhist, ms)
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
                solution_alleles = set(a for a, p in proportional_solution_alleles if p > decimal.Decimal('0.3'))  # small prop cases will be disregarded in non-proportional calling
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
def test_prf_bi_proportional_calling(histograms_and_calling_solutions_d, prf_simcorbipropschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = prf_simcorbipropschema.call_ms_hist(dbhist, ms)
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


@pytest.mark.django_db
def test_prf_bi_proportional_highest_peak_calling(histograms_and_calling_solutions_d, prf_simcorbiprophighpeakschema):
    for amp_id, ms_dict in histograms_and_calling_solutions_d.items():
        for ms_id, histograms_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for proportional_solution_alleles, dbhist in histograms_dict.items():
                bcca = prf_simcorbiprophighpeakschema.call_ms_hist(dbhist, ms)
                assert isinstance(bcca, BestCorrelationProportionalHighestPeakCalledAlleles)
                bcca = BestCorrelationProportionalHighestPeakCalledAlleles.objects.get(pk=bcca.pk)  # get rid of possible overriden field
                pmas = bcca.genotypes.proportionalmicrosatellitealleleset
                assert isinstance(pmas, ProportionalMicrosatelliteAlleleSet)
                proportional_alleles = pmas.alleles
                result_allels = set(a for a, p in proportional_alleles)
                solution_alleles = set(a for a, p in proportional_solution_alleles)
                assert result_allels == solution_alleles  # assert called alleles against solution alleles
                allowed_genotypes = set(frozenset([(a, p) for a, p in aap if p > 0]) for aap, c in prf_simcorbiprophighpeakschema.alleles_and_cycles)
                if proportional_solution_alleles in allowed_genotypes:
                    assert proportional_alleles == proportional_solution_alleles  # calling should be identical to solution
                else:
                    small_prop_allele, large_prop_allele = sorted(proportional_alleles, key=lambda t: t[1])
                    small_prop_solution_allele, large_prop_solution_allele = sorted(proportional_solution_alleles, key=lambda t: t[1])
                    allele, prop = small_prop_solution_allele
                    small_prop_solution_allele = allele, prop + decimal.Decimal('0.1')
                    assert small_prop_allele == small_prop_solution_allele
                    allele, prop = large_prop_solution_allele
                    large_prop_solution_allele = allele, prop - decimal.Decimal('0.1')
                    assert large_prop_allele == large_prop_solution_allele
