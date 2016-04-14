import pytest

from targeted_enrichment.amplicons.models import PlainTargetedAmplicon

from tests.targeted_enrichment.planning.conftest import *


@pytest.fixture()
def pu_28727(ugs_28727_left, ugs_28727_right):
    pu = PlainTargetedAmplicon(
        id=28727,
        left_ugs = ugs_28727_left,
        right_ugs = ugs_28727_right,
    )
    pu.infer_slice()  # FIXME
    pu.save()
    return pu


@pytest.fixture()
def pu_28734(ugs_28734_left, ugs_28734_right):
    pu = PlainTargetedAmplicon(
        id=28734,
        left_ugs = ugs_28734_left,
        right_ugs = ugs_28734_right,
    )
    pu.infer_slice()  # FIXME
    pu.save()
    return pu


@pytest.fixture()
def requires_amplicons(pu_28727, pu_28734):
    pass
