import pytest


@pytest.mark.django_db
def test_assembly(hg19_assembly):
    assert hg19_assembly.friendly_name == "hg19"


@pytest.mark.django_db
def test_chromosome(hg19_chromosome):
    assert hg19_chromosome.name == "X"
