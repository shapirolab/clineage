import pytest
from models import PCR1Multiplex, Panel, PCR1MultiplexCollection


@pytest.fixture()
def pcr1multiplex(db):
    pcr1m = PCR1Multiplex.objects.create(
        name='test PCR1Multiplex'
    )
    return pcr1m

@pytest.fixture()
def panel(db):
    p = Panel.objects.create(
        name='test Panel'
    )
    return p


@pytest.fixture()
def pcr1multiplexcollection(panel):
    pcr1mc = PCR1MultiplexCollection.objects.create(
        panel=panel
    )
    return pcr1mc


@pytest.mark.django_db
def test_pcr1multiplex(pcr1multiplex):
    assert pcr1multiplex.name == 'test PCR1Multiplex'


@pytest.mark.django_db
def test_pcr1multiplexcollection(pcr1multiplexcollection):
    assert pcr1multiplexcollection.panel.name == 'test Panel'
