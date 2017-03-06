import pytest

from targeted_enrichment.amplicons.models import PlainTargetedAmplicon, AmpliconCollection, UMITargetedAmplicon

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
def pu_snp_example_1(ugs_snp_1_left, ugs_snp_1_right, slice_snp_1_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_snp_1_amplicon,
        id=5,
        left_ugs=ugs_snp_1_left,
        right_ugs=ugs_snp_1_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def pu_snp_example_2(ugs_snp_2_left, ugs_snp_2_right, slice_snp_2_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_snp_2_amplicon,
        id=6,
        left_ugs=ugs_snp_2_left,
        right_ugs=ugs_snp_2_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def pu_snp_example_3(ugs_snp_3_left, ugs_snp_3_right, slice_snp_3_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_snp_3_amplicon,
        id=8,
        left_ugs=ugs_snp_3_left,
        right_ugs=ugs_snp_3_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def pu_snp_example_4(ugs_snp_4_left, ugs_snp_4_right, slice_snp_4_amplicon):
    pu = PlainTargetedAmplicon.objects.create(
        slice=slice_snp_4_amplicon,
        id=31584,
        left_ugs=ugs_snp_4_left,
        right_ugs=ugs_snp_4_right,
    )
    # So our objects don't have "special" objects in fields
    pu = PlainTargetedAmplicon.objects.get(pk=pu.pk)
    return pu


@pytest.fixture()
def utm_snp_example_1(ugs_snp_1_left, ugs_snp_1_right, slice_snp_1_amplicon):
    pu = UMITargetedAmplicon.objects.create(
        slice=slice_snp_1_amplicon,
        id=7,
        left_ugs=ugs_snp_1_left,
        right_ugs=ugs_snp_1_right,
        umi_length=3,
    )
    # So our objects don't have "special" objects in fields
    pu = UMITargetedAmplicon.objects.get(pk=pu.pk)
    return pu

@pytest.fixture()
def requires_amplicons(pu_28727, pu_28734, pu_adj_ms_1, pu_adj_ms_2):
    pass


@pytest.fixture()
def requires_snp_amplicons(pu_snp_example_1, pu_snp_example_2, pu_snp_example_3, pu_snp_example_4, utm_snp_example_1):
    pass

@pytest.fixture()
def requires_snp_targets(slice_snp_1_target_a, slice_snp_1_target_b, slice_snp_2_target_a):
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
    ac.amplicons = list(PlainTargetedAmplicon.objects.filter(id__in=[1, 3, 4]))
    ac.save()
    # So our objects don't have "special" objects in fields
    ac = AmpliconCollection.objects.get(pk=ac.pk)
    return ac


@pytest.fixture()
def snp_1234(requires_snp_amplicons):
    snp_example = AmpliconCollection.objects.create()
    l = list()
    l.append(PlainTargetedAmplicon.objects.get(id=5))
    l.append(PlainTargetedAmplicon.objects.get(id=6))
    l.append(PlainTargetedAmplicon.objects.get(id=31584))
    l.append(PlainTargetedAmplicon.objects.get(id=8))
    l.append(UMITargetedAmplicon.objects.get(id=7))
    snp_example.amplicons = l
    snp_example.save()
    # So our objects don't have "special" objects in fields
    ac = AmpliconCollection.objects.get(pk=snp_example.pk)
    return snp_example
