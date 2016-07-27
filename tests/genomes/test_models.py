import pytest


@pytest.mark.django_db
def test_assembly(hg19_assembly):
    assert str(hg19_assembly) == "hg19"


@pytest.mark.django_db
def test_chromosome(hg19_chromosome):
    assert str(hg19_chromosome) == "hg19:X"


@pytest.mark.django_db
def test_dnaslice(slice_28727_left):
    assert str(slice_28727_left) == "X:81316094-81316116"


@pytest.mark.skipif(pytest.config.getoption("nomigrations"), reason="No migrations, no view.")
@pytest.mark.django_db
def test_dnasliceintersection(slice_28727_target_a, slice_28727_target_b, slice_28727_amplicon):
    assert set(slice_28727_amplicon.contains.all()) == set([slice_28727_target_a, slice_28727_target_b])
    assert set(slice_28727_amplicon.contained.all()) == set()
    assert set(slice_28727_target_a.contains.all()) == set()
    assert set(slice_28727_target_a.contained.all()) == set([slice_28727_amplicon])
    assert set(slice_28727_target_b.contains.all()) == set()
    assert set(slice_28727_target_b.contained.all()) == set([slice_28727_amplicon])
