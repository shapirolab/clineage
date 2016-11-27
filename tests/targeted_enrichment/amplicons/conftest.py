import pytest

from targeted_enrichment.amplicons.models import PlainTargetedAmplicon, AmpliconCollection

from tests.targeted_enrichment.planning.conftest import *
from tests.genomes.conftest import *


@pytest.fixture()
def pu_28727(ugs_28727_left, ugs_28727_right, slice_28727_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_28727_amplicon,
        id=1,
        left_ugs=ugs_28727_left,
        right_ugs=ugs_28727_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def pu_28734(ugs_28734_left, ugs_28734_right, slice_28734_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_28734_amplicon,
        id=2,
        left_ugs=ugs_28734_left,
        right_ugs=ugs_28734_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def pu_adj_ms_1(ugs_adj_ms_1_left, ugs_adj_ms_1_right, slice_adj_ms_1_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_adj_ms_1_amplicon,
        id=3,
        left_ugs=ugs_adj_ms_1_left,
        right_ugs=ugs_adj_ms_1_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def pu_adj_ms_2(ugs_adj_ms_2_left, ugs_adj_ms_2_right, slice_adj_ms_2_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_adj_ms_2_amplicon,
        id=4,
        left_ugs=ugs_adj_ms_2_left,
        right_ugs=ugs_adj_ms_2_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def requires_amplicons(pu_28727, pu_28734, pu_adj_ms_1, pu_adj_ms_2):
    pass


@pytest.fixture()
def amplicon_collection(requires_amplicons):
    ac = AmpliconCollection.objects.create()
    ac.amplicons = list(PlainTargetedAmplicon.objects.all())
    ac.save()
    # So our objects don't have "special" objects in fields
    ac = AmpliconCollection.objects.get(pk=ac.pk)
    return ac


@pytest.fixture()
def ac_134(requires_amplicons):
    ac = AmpliconCollection.objects.create()
    ac.amplicons = list(PlainTargetedAmplicon.objects.filter(id__in=[1,3,4]))
    ac.save()
    # So our objects don't have "special" objects in fields
    ac = AmpliconCollection.objects.get(pk=ac.pk)
    return ac