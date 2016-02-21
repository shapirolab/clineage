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


@pytest.mark.django_db
def test_taxa(hg19_assembly):
    assert human_taxa.friendly_name == "Human"