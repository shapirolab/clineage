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
