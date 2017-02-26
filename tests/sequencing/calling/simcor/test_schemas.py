import pytest

@pytest.mark.django_db
def test_mono_schema(minimalsimcormonoschema):
    assert set(minimalsimcormonoschema.alleles_and_cycles) == set([
        (15, 20), (16, 20), (15, 21), (16, 21),
        (15, 22), (16, 22), (15, 23), (16, 23)
    ])


@pytest.mark.django_db
def test_bi_schema(minimalsimcorbischema):
    assert minimalsimcorbischema.allele_number == 2
    assert set(minimalsimcorbischema.alleles_and_cycles) == set([
        (frozenset([15, 16]), 20), (frozenset([15, ]), 20), (frozenset([16, ]), 20),
        (frozenset([15, 16]), 21), (frozenset([15, ]), 21), (frozenset([16, ]), 21),
        (frozenset([15, 16]), 22), (frozenset([15, ]), 22), (frozenset([16, ]), 22),
        (frozenset([15, 16]), 23), (frozenset([15, ]), 23), (frozenset([16, ]), 23),
    ])


@pytest.mark.django_db
def test_proportional_bi_schema(minimalsimcorbipropschema):
    assert minimalsimcorbipropschema.allele_number == 2
    assert set(minimalsimcorbipropschema.alleles_and_cycles) == set([
        (frozenset([(15, 16]), 20), (frozenset([15, ]), 20), (frozenset([16, ]), 20),
    ])