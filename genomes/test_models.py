import pytest
from models import Assembly, Chromosome
from misc.test_models import human_taxa

@pytest.fixture()
def hg19_assembly(human_taxa):
    a = Assembly.objects.create(
        taxa=human_taxa,
        name='February 2009 Homo sapiens (GRCh37)',
        friendly_name='hg19',
    )
    return a


@pytest.fixture()
def hg19_chromosome(hg19_assembly):
    c = Chromosome.objects.create(
        assembly=hg19_assembly,
        name='X',
        sequence_length=155270560,
        cyclic=False,
    )
    return c


@pytest.mark.django_db
def test_assembly(hg19_assembly):
    assert hg19_assembly.friendly_name == "hg19"


@pytest.mark.django_db
def test_chromosome(hg19_chromosome):
    assert hg19_chromosome.name == "X"