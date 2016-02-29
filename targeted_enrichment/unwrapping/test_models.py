import pytest
from models import PCR1Unwrapper

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

from ..reagents.test_models import ter_28727, ter_28734


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




@pytest.mark.django_db
def test_pu_28727(pu_28727, ter_28727):
    assert pu_28727.ter == ter_28727
    assert pu_28727.left_margin.seq == "ATTTACTATGCCATGCTGCTGCT"
