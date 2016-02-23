import pytest
from models import PCR1PrimerPairTER

from misc.test_models import human_taxa
from genomes.test_models import hg19_assembly, hg19_chromosome, \
    slice_28727_left, slice_28727_right, slice_28727_target_a, slice_28727_target_b,\
    slice_28734_left, slice_28734_right, slice_28734_target_a
from ..planning.test_models import te_28727, te_28734, \
    ugs_28727_left, ugs_28727_right, ugs_28734_left, ugs_28734_right, \
    ms_28727_a, ms_28727_b, ms_28734_a
from primers.synthesis.test_models import primer_28727_left, primer_28727_right, \
    primer_28734_left, primer_28734_right
from primers.parts.test_models import illuminareadingadaptor1cuts, illuminareadingadaptor2cuts, \
    illuminareadingadaptor1, illuminareadingadaptor2


@pytest.fixture()
def ter_28727(te_28727, primer_28727_left, primer_28727_right):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_28727,
        passed_validation=False,
        left_primer=primer_28727_left,
        right_primer=primer_28727_right,
    )
    return ter


@pytest.fixture()
def ter_28734(te_28734, primer_28734_left, primer_28734_right):
    ter = PCR1PrimerPairTER.objects.create(
        te=te_28734,
        passed_validation=True,
        left_primer=primer_28734_left,
        right_primer=primer_28734_right,
    )
    return ter


@pytest.mark.django_db
def test_ter(ter_28727):
    assert ter_28727.left_primer.sequence.seq == "CTACACGACGCTCTTCCGATCTATTTACTATGCCATGCTGCTGCT"
    assert ter_28727.te.targets.get(name='X_81316201_81316236').slice.start_pos == 81316201