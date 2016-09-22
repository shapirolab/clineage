import pytest

from targeted_enrichment.planning.models import UGSPlus, UGSMinus, \
    TargetEnrichment, Microsatellite, RestrictionEnzyme

from tests.genomes.conftest import *


@pytest.fixture()
def ugs_28727_left(slice_28727_left):
    ugsp = UGSPlus.objects.create(
        slice=slice_28727_left,
    )
    # So our objects don't have "special" objects in fields
    ugsp = UGSPlus.objects.get(pk=ugsp.pk)
    return ugsp


@pytest.fixture()
def ugs_28727_right(slice_28727_right):
    ugsm = UGSMinus.objects.create(
        slice=slice_28727_right,
    )
    # So our objects don't have "special" objects in fields
    ugsm = UGSMinus.objects.get(pk=ugsm.pk)
    return ugsm


@pytest.fixture()
def ms_28727_a(slice_28727_target_a):
    ms = Microsatellite.objects.create(
        id=1,
        name='X_81316131_81316199',
        slice=slice_28727_target_a,
        repeat_unit_len=3,
        repeat_unit_type='AAG',
        repeat_number=23,
        repeat_unit_ref_seq='TCT',
        planning_version=0,
    )
    # So our objects don't have "special" objects in fields
    ms = Microsatellite.objects.get(pk=ms.pk)
    return ms


@pytest.fixture()
def ms_28727_b(slice_28727_target_b):
    ms = Microsatellite.objects.create(
        id=2,
        name='X_81316201_81316236',
        slice=slice_28727_target_b,
        repeat_unit_len=3,
        repeat_unit_type='AGC',
        repeat_number=12,
        repeat_unit_ref_seq='CTG',
        planning_version=0,
    )
    # So our objects don't have "special" objects in fields
    ms = Microsatellite.objects.get(pk=ms.pk)
    return ms


@pytest.fixture()
def te_28727(hg19_chromosome, ugs_28727_left, ugs_28727_right):
    te = TargetEnrichment.objects.create(
        chromosome=hg19_chromosome,
        left=ugs_28727_left,
        right=ugs_28727_right,
        planning_version=1,
    )
    # So our objects don't have "special" objects in fields
    te = TargetEnrichment.objects.get(pk=te.pk)
    return te


@pytest.fixture()
def ugs_28734_left(slice_28734_left):
    ugsp = UGSPlus.objects.create(
        slice=slice_28734_left,
    )
    # So our objects don't have "special" objects in fields
    ugsp = UGSPlus.objects.get(pk=ugsp.pk)
    return ugsp


@pytest.fixture()
def ugs_28734_right(slice_28734_right):
    ugsm = UGSMinus.objects.create(
        slice=slice_28734_right,
    )
    # So our objects don't have "special" objects in fields
    ugsm = UGSMinus.objects.get(pk=ugsm.pk)
    return ugsm


@pytest.fixture()
def ugs_adj_ms_1_left(slice_adj_ms_1_left):
    ugsp = UGSPlus.objects.create(
        slice=slice_adj_ms_1_left,
    )
    # So our objects don't have "special" objects in fields
    ugsp = UGSPlus.objects.get(pk=ugsp.pk)
    return ugsp


@pytest.fixture()
def ugs_adj_ms_1_right(slice_adj_ms_1_right):
    ugsm = UGSMinus.objects.create(
        slice=slice_adj_ms_1_right,
    )
    # So our objects don't have "special" objects in fields
    ugsm = UGSMinus.objects.get(pk=ugsm.pk)
    return ugsm


@pytest.fixture()
def ugs_adj_ms_2_left(slice_adj_ms_2_left):
    ugsp = UGSPlus.objects.create(
        slice=slice_adj_ms_2_left,
    )
    # So our objects don't have "special" objects in fields
    ugsp = UGSPlus.objects.get(pk=ugsp.pk)
    return ugsp


@pytest.fixture()
def ugs_adj_ms_2_right(slice_adj_ms_2_right):
    ugsm = UGSMinus.objects.create(
        slice=slice_adj_ms_2_right,
    )
    # So our objects don't have "special" objects in fields
    ugsm = UGSMinus.objects.get(pk=ugsm.pk)
    return ugsm


@pytest.fixture()
def ms_28734_a(slice_28734_target_a):
    ms = Microsatellite.objects.create(
        id=3,
        name='X_54384788_54384805',
        slice=slice_28734_target_a,
        repeat_unit_len=3,
        repeat_unit_type='AAG',
        repeat_number=6,
        repeat_unit_ref_seq='AGA',
        planning_version=0,
    )
    # So our objects don't have "special" objects in fields
    ms = Microsatellite.objects.get(pk=ms.pk)
    return ms


@pytest.fixture()
def te_28734(hg19_chromosome, ugs_28734_left, ugs_28734_right):
    te = TargetEnrichment.objects.create(
        chromosome=hg19_chromosome,
        left=ugs_28734_left,
        right=ugs_28734_right,
        planning_version=1,
    )
    # So our objects don't have "special" objects in fields
    te = TargetEnrichment.objects.get(pk=te.pk)
    return te


@pytest.fixture()
def ttaa_restriction_site_sample():
    r_site = RestrictionEnzyme.objects.create(
        name='MseI',
        sequence='TTAA',
        cut_delta=2,
        sticky_bases=2,
        sequence_len=4,
    )
    # So our objects don't have "special" objects in fields
    r_site = RestrictionEnzyme.objects.get(name=r_site.name)
    return r_site


@pytest.fixture()
def requires_microsatellites(ms_28727_a, ms_28727_b, ms_28734_a, ms_adj_ms_1_a, ms_adj_ms_2_a, ms_adj_ms_2_b):
    pass

