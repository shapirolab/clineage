import pytest
from models import PCR1PlusPrimer, PCR1MinusPrimer

from misc.test_models import human_taxa
from genomes.test_models import hg19_assembly, hg19_chromosome, \
    slice_28727_left, slice_28727_right, \
    slice_28734_left, slice_28734_right
from targeted_enrichment.planning.test_models import ugs_28727_left, ugs_28734_right, ugs_28734_left, ugs_28727_right
from ..parts.test_models import illuminareadingadaptor1cuts, illuminareadingadaptor2cuts, \
    illuminareadingadaptor1, illuminareadingadaptor2


@pytest.fixture()
def primer_28727_left(ugs_28727_left, illuminareadingadaptor1cuts):
    p = PCR1PlusPrimer.objects.create(
        name='chrX_81316201_81316236_ACG_fwd',
        ugs=ugs_28727_left,
        irac=illuminareadingadaptor1cuts,
    )
    return p


@pytest.fixture()
def primer_28727_right(ugs_28727_right, illuminareadingadaptor2cuts):
    p = PCR1MinusPrimer.objects.create(
        name='Seq19272_rev',
        ugs=ugs_28727_right,
        irac=illuminareadingadaptor2cuts,
    )
    return p


@pytest.fixture()
def primer_28734_left(ugs_28734_left, illuminareadingadaptor1cuts):
    p = PCR1PlusPrimer.objects.create(
        name='Seq20808_fwd',
        ugs=ugs_28734_left,
        irac=illuminareadingadaptor1cuts,
    )
    return p


@pytest.fixture()
def primer_28734_right(ugs_28734_right, illuminareadingadaptor2cuts):
    p = PCR1MinusPrimer.objects.create(
        name='Seq20808_rev',
        ugs=ugs_28734_right,
        irac=illuminareadingadaptor2cuts,
    )
    return p



@pytest.mark.django_db
def test_primer_28727_left_head(primer_28727_left):
    assert primer_28727_left.head.seq == "ATTTACTATGCCATGCTGCTGCT"


@pytest.mark.django_db
def test_primer_28727_right_head(primer_28727_right):
    assert primer_28727_right.head.seq == "GGCATCTGTTCTTGTTTGCACAT"


@pytest.mark.django_db
def test_primer_28727_left_tail(primer_28727_left):
    assert primer_28727_left.tail.seq == "CTACACGACGCTCTTCCGATCT"


@pytest.mark.django_db
def test_primer_28727_right_tail(primer_28727_right):
    assert primer_28727_right.tail.seq == "CAGACGTGTGCTCTTCCGATCT"


@pytest.mark.django_db
def test_primer_28727_left(primer_28727_left):
    assert primer_28727_left.sequence.seq == "CTACACGACGCTCTTCCGATCTATTTACTATGCCATGCTGCTGCT"


@pytest.mark.django_db
def test_primer_28727_right(primer_28727_right):
    assert primer_28727_right.sequence.seq == "CAGACGTGTGCTCTTCCGATCTGGCATCTGTTCTTGTTTGCACAT"


@pytest.mark.django_db
def test_primer_28734_left(primer_28734_left):
    assert primer_28734_left.sequence.seq == "CTACACGACGCTCTTCCGATCTAAAGGCTTCTCCCCATTCCAAAG"


@pytest.mark.django_db
def test_primer_28734_right(primer_28734_right):
    assert primer_28734_right.sequence.seq == "CAGACGTGTGCTCTTCCGATCTGGAAGTAGTGTGTGCTTGGACTT"