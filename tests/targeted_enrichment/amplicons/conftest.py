import pytest

from targeted_enrichment.amplicons.models import PCR1Amplicon

from tests.targeted_enrichment.reagents.conftest import *


@pytest.fixture()
def pu_28727(ter_28727):
    pu = PCR1Amplicon(
        id=28727,
        ter=ter_28727,
    )
    pu.infer_slice()
    pu.save()
    return pu


@pytest.fixture()
def pu_28734(ter_28734):
    pu = PCR1Amplicon(
        id=28734,
        ter=ter_28734,
    )
    pu.infer_slice()
    pu.save()
    return pu


@pytest.fixture()
def requires_amplicons(pu_28727, pu_28734):
    pass
