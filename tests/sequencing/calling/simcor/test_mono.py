import pytest
from sequencing.calling.simcor.calling import call_microsatellite_histogram


@pytest.mark.django_db
def test_mono_calling(histograms_and_calling_solutions, simcorschema):
    for solution_alleles, dbhist, ms in histograms_and_calling_solutions:
        # TODO: make this a function of calling schemas
        bcca = call_microsatellite_histogram(dbhist, ms, simcorschema)
        alleles = bcca.genotypes.alleles
        assert alleles == set(a for p, a in solution_alleles)  # assert called alleles against solution alleles
