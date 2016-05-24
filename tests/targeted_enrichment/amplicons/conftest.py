import pytest

from targeted_enrichment.amplicons.models import PlainTargetedAmplicon

from tests.targeted_enrichment.planning.conftest import *
from tests.genomes.conftest import *


@pytest.fixture()
def pu_28727(ugs_28727_left, ugs_28727_right, slice_28727_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_28727_amplicon,
        id=28727,
        left_ugs = ugs_28727_left,
        right_ugs = ugs_28727_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def pu_28734(ugs_28734_left, ugs_28734_right, slice_28734_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_28734_amplicon,
        id=28734,
        left_ugs = ugs_28734_left,
        right_ugs = ugs_28734_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def requires_amplicons(pu_28727, pu_28734):
    pass
