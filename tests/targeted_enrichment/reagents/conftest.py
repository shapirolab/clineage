import pytest

from targeted_enrichment.reagents.models import PCR1PrimerPairTER

from tests.targeted_enrichment.planning.conftest import *
from tests.targeted_enrichment.amplicons.conftest import *
from tests.primers.synthesis.conftest import *


@pytest.fixture()
def ter_28727(te_28727, primer_28727_left, primer_28727_right, pu_28727):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_28727,
        passed_validation=False,
        left_primer=primer_28727_left,
        right_primer=primer_28727_right,
        amplicon=pu_28727,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def ter_28734(te_28734, primer_28734_left, primer_28734_right, pu_28734):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_28734,
        passed_validation=True,
        left_primer=primer_28734_left,
        right_primer=primer_28734_right,
        amplicon=pu_28734,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def ter_adj_ms_1(te_adj_ms_1, primer_adj_ms_1_left, primer_adj_ms_1_right, pu_adj_ms_1):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_adj_ms_1,
        passed_validation=True,
        left_primer=primer_adj_ms_1_left,
        right_primer=primer_adj_ms_1_right,
        amplicon=pu_adj_ms_1,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def ter_adj_ms_2(te_adj_ms_2, primer_adj_ms_2_left, primer_adj_ms_2_right, pu_adj_ms_2):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_adj_ms_2,
        passed_validation=True,
        left_primer=primer_adj_ms_2_left,
        right_primer=primer_adj_ms_2_right,
        amplicon=pu_adj_ms_2,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def ter_snp_1(te_snp_1, primer_snp_1_left, primer_snp_1_right, pu_snp_example_1):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_snp_1,
        passed_validation=False,
        left_primer=primer_snp_1_left,
        right_primer=primer_snp_1_right,
        amplicon=pu_snp_example_1,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def ter_snp_2(te_snp_2, primer_snp_2_left, primer_snp_2_right, pu_snp_example_2):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_snp_2,
        passed_validation=False,
        left_primer=primer_snp_2_left,
        right_primer=primer_snp_2_right,
        amplicon=pu_snp_example_2,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def ter_snp_3(te_snp_3, primer_snp_3_left, primer_snp_3_right, pu_snp_example_3):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_snp_3,
        passed_validation=False,
        left_primer=primer_snp_3_left,
        right_primer=primer_snp_3_right,
        amplicon=pu_snp_example_3,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def ter_snp_4(te_snp_4, primer_snp_4_left, primer_snp_4_right, pu_snp_example_4):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_snp_4,
        passed_validation=False,
        left_primer=primer_snp_4_left,
        right_primer=primer_snp_4_right,
        amplicon=pu_snp_example_4,
    )
    # So our objects don't have "special" objects in fields
    ter = PCR1PrimerPairTER.objects.get(pk=ter.pk)
    return ter


@pytest.fixture()
def all_ters(ter_28727, ter_28734, ter_adj_ms_1, ter_adj_ms_2, ter_snp_1, ter_snp_2, ter_snp_3, ter_snp_4):
    return [ter_28727, ter_28734, ter_adj_ms_1, ter_adj_ms_2, ter_snp_1, ter_snp_2, ter_snp_3, ter_snp_4]