import pytest


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
