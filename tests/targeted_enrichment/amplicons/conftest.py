import pytest

from targeted_enrichment.amplicons.models import PCR1Amplicon

from tests.targeted_enrichment.reagents.conftest import *


@pytest.fixture()
def pu_28727(ter_28727):
    pu = PCR1Amplicon(
        id=1,
        ter=ter_28727,
    )
    pu.infer_slice()
    pu.save()
    return pu


@pytest.fixture()
def pu_28734(ter_28734):
    pu = PCR1Amplicon(
        id=2,
        ter=ter_28734,
    )
    pu.infer_slice()
    pu.save()
    return pu


@pytest.fixture()
def amplicon_d(pu_28727, pu_28734):
    return {
        "28727": pu_28727,
        "28734": pu_28734,
    }


@pytest.fixture()
def amplicon_d_r(amplicon_d):
    return {v: k for k, v in amplicon_d.iteritems()}
