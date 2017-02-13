import pytest

from lib_prep.multiplexes.models import PCR1Multiplex, PCR1Panel, OM6Panel, OM6Oligomix

from tests.targeted_enrichment.planning.conftest import *
from tests.targeted_enrichment.reagents.conftest import *


@pytest.fixture()
def pcr1multiplex(all_ters):
    pcr1m = PCR1Multiplex.objects.create(
        name='test PCR1Multiplex'
    )
    pcr1m.ters = all_ters
    # So our objects don't have "special" objects in fields
    pcr1m = PCR1Multiplex.objects.get(pk=pcr1m.pk)
    return pcr1m


@pytest.fixture()
def pcr1panel(pcr1multiplex, amplicon_collection):
    pcr1mc = PCR1Panel.objects.create(
        name='test PCR1Panel',
        amplicon_collection=amplicon_collection,
    )
    pcr1mc.mpxs = [pcr1multiplex]
    # So our objects don't have "special" objects in fields
    pcr1mc = PCR1Panel.objects.get(pk=pcr1mc.pk)
    return pcr1mc

@pytest.fixture()
def padlockmix(all_ters_padlock):
    pdmix = OM6Oligomix.objects.create(
        name='test padlockmix'
    )
    pdmix.ters = all_ters_padlock
    # So our objects don't have "special" objects in fields
    pdmix = OM6Oligomix.objects.get(pk=pdmix.pk)
    return pdmix


@pytest.fixture()
def padlockpanel(padlockmix, amplicon_collection):
    pdpanel = OM6Panel.objects.create(
        name='test_padlockpanel',
        amplicon_collection=amplicon_collection,
    )
    pdpanel.mixs = [padlockmix]
    # So our objects don't have "special" objects in fields
    pdpanel = OM6Panel.objects.get(pk=pdpanel.pk)
    return pdpanel
