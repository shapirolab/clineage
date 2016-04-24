import pytest

from lib_prep.multiplexes.models import PCR1Multiplex, PCR1Panel

from tests.targeted_enrichment.planning.conftest import *
from tests.targeted_enrichment.reagents.conftest import *


@pytest.fixture()
def pcr1multiplex(ter_28727, ter_28734):
    pcr1m = PCR1Multiplex.objects.create(
        name='test PCR1Multiplex'
    )
    pcr1m.ters = [ter_28727, ter_28734]
    return pcr1m


@pytest.fixture()
def pcr1panel(pcr1multiplex):
    pcr1mc = PCR1Panel.objects.create(
        name='test PCR1Panel'
    )
    pcr1mc.mpxs = [pcr1multiplex]
    return pcr1mc
