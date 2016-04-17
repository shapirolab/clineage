import pytest

from targeted_enrichment.planning.models import UGSPlus, UGSMinus, \
    TargetEnrichment, Microsatellite, PhasedMicrosatellites

from tests.genomes.conftest import *


@pytest.fixture()
def ugs_28727_left(slice_28727_left):
    ugsp = UGSPlus.objects.create(
        slice=slice_28727_left,
    )
    return ugsp


@pytest.fixture()
def ugs_28727_right(slice_28727_right):
    ugsm = UGSMinus.objects.create(
        slice=slice_28727_right,
    )
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
    )
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
    )
    return ms


@pytest.fixture()
def te_28727(hg19_chromosome, ugs_28727_left, ugs_28727_right):
    te = TargetEnrichment.objects.create(
        chromosome=hg19_chromosome,
        left=ugs_28727_left,
        right=ugs_28727_right,
        planning_version=1,
    )
    return te


@pytest.fixture()
def pms_28727(ms_28727_a, ms_28727_b, slice_28727_amplicon):
    pms = PhasedMicrosatellites.objects.create(
        slice=slice_28727_amplicon,
        planning_version=0,
    )
    pms.microsatellites = [ms_28727_a, ms_28727_b]
    return pms


@pytest.fixture()
def ugs_28734_left(slice_28734_left):
    ugsp = UGSPlus.objects.create(
        slice=slice_28734_left,
    )
    return ugsp


@pytest.fixture()
def ugs_28734_right(slice_28734_right):
    ugsm = UGSMinus.objects.create(
        slice=slice_28734_right,
    )
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
    )
    return ms


@pytest.fixture()
def te_28734(hg19_chromosome, ugs_28734_left, ugs_28734_right):
    te = TargetEnrichment.objects.create(
        chromosome=hg19_chromosome,
        left=ugs_28734_left,
        right=ugs_28734_right,
        planning_version=1,
    )
    return te


@pytest.fixture()
def pms_28734(ms_28734_a, slice_28734_amplicon):
    pms = PhasedMicrosatellites.objects.create(
        slice=slice_28734_amplicon,
        planning_version=0,
    )
    pms.microsatellites = [ms_28734_a]
    return pms


@pytest.fixture()
def requires_pmss(pms_28727, pms_28734):
    pass
