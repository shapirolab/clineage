import pytest

from targeted_enrichment.unwrapping.models import PCR1Unwrapper

from tests.targeted_enrichment.reagents.conftest import *


@pytest.fixture()
def pu_28727(ter_28727):
    pu = PCR1Unwrapper(
        ter=ter_28727
    )
    pu.infer_slice()
    pu.save()
    return pu


@pytest.fixture()
def pu_28734(ter_28734):
    pu = PCR1Unwrapper(
        ter=ter_28734
    )
    pu.infer_slice()
    pu.save()
    return pu

@pytest.fixture()
def require_unwrappers(pu_28727, pu_28734):
    pass