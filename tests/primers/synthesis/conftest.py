import pytest

from primers.synthesis.models import PCR1PlusPrimer, PCR1MinusPrimer

from tests.targeted_enrichment.planning.conftest import *
from tests.primers.parts.conftest import *


@pytest.fixture()
def primer_28727_left(ugs_28727_left, illuminareadingadaptor1cuts):
    p = PCR1PlusPrimer.objects.create(
        name='chrX_81316201_81316236_ACG_fwd',
        ugs=ugs_28727_left,
        irac=illuminareadingadaptor1cuts,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_28727_right(ugs_28727_right, illuminareadingadaptor2cuts):
    p = PCR1MinusPrimer.objects.create(
        name='Seq19272_rev',
        ugs=ugs_28727_right,
        irac=illuminareadingadaptor2cuts,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_28734_left(ugs_28734_left, illuminareadingadaptor1cuts):
    p = PCR1PlusPrimer.objects.create(
        name='Seq20808_fwd',
        ugs=ugs_28734_left,
        irac=illuminareadingadaptor1cuts,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_28734_right(ugs_28734_right, illuminareadingadaptor2cuts):
    p = PCR1MinusPrimer.objects.create(
        name='Seq20808_rev',
        ugs=ugs_28734_right,
        irac=illuminareadingadaptor2cuts,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p

