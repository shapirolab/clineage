import pytest

from primers.synthesis.models import PCR1PlusPrimer, PCR1MinusPrimer, OM6Padlock
from primers.synthesis.models import OM6Padlock, OM6Prep

from tests.targeted_enrichment.planning.conftest import *
from tests.primers.parts.conftest import *
from tests.amplicons_handling.conftest import *


@pytest.fixture()
def primer_28727_left(ugs_28727_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='chrX_81316201_81316236_ACG_fwd',
        ugs=ugs_28727_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_28727_right(ugs_28727_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='Seq19272_rev',
        ugs=ugs_28727_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_28734_left(ugs_28734_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='Seq20808_fwd',
        ugs=ugs_28734_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_28734_right(ugs_28734_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='Seq20808_rev',
        ugs=ugs_28734_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_adj_ms_1_left(ugs_adj_ms_1_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='Seq00044_fwd',
        ugs=ugs_adj_ms_1_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_adj_ms_1_right(ugs_adj_ms_1_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='Seq00044_rev',
        ugs=ugs_adj_ms_1_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_adj_ms_2_left(ugs_adj_ms_2_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='Seq01291_fwd',
        ugs=ugs_adj_ms_2_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_adj_ms_2_right(ugs_adj_ms_2_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='Seq01291_rev',
        ugs=ugs_adj_ms_2_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_1_left(ugs_snp_1_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='chr1_115256468_115256612_fwd',
        ugs=ugs_snp_1_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_1_right(ugs_snp_1_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='chr1_115256468_115256612_rev',
        ugs=ugs_snp_1_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_2_left(ugs_snp_2_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='chr1_115256468_115256487_fwd',
        ugs=ugs_snp_2_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_2_right(ugs_snp_2_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='chr10_10341355_10341542_rev',
        ugs=ugs_snp_2_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_3_left(ugs_snp_3_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='chr4_106155054_106155263_fwd',
        ugs=ugs_snp_3_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_3_right(ugs_snp_3_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='chr4_106155054_106155263_rev',
        ugs=ugs_snp_3_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_4_left(ugs_snp_4_left, illuminareadingadaptor1fortail):
    p = PCR1PlusPrimer.objects.create(
        name='chr21_36259181_36259442_fwd',
        ugs=ugs_snp_4_left,
        iraft=illuminareadingadaptor1fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1PlusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def primer_snp_4_right(ugs_snp_4_right, illuminareadingadaptor2fortail):
    p = PCR1MinusPrimer.objects.create(
        name='chr21_36259181_36259442_rev',
        ugs=ugs_snp_4_right,
        iraft=illuminareadingadaptor2fortail,
    )
    # So our objects don't have "special" objects in fields
    p = PCR1MinusPrimer.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def padlock_snp_1(ugs_snp_4_left,ugs_snp_4_right,illuminareadingadaptor1fortail,illuminareadingadaptor2fortail):
    p = OM6Padlock.objects.create(
        left_ugs=ugs_snp_4_left,
        right_ugs=ugs_snp_4_right,
        ira1ft=illuminareadingadaptor1fortail,
        ira2ft=illuminareadingadaptor2fortail,
        umi_length=3,
    )
    # So our objects don't have "special" objects in fields
    p = OM6Padlock.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def padlock_prep(padlock_snp_1, padlock_prep_primers):
    om6p = OM6Prep.objects.create(
        padlock=padlock_snp_1,
        primers=padlock_prep_primers
    )
    # So our objects don't have "special" objects in fields
    om6p = OM6Prep.objects.get(pk=om6p.pk)
    return om6p
