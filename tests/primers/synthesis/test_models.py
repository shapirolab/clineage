import pytest


@pytest.mark.django_db
def test_primer_28727_left_head(primer_28727_left):
    assert primer_28727_left.head.seq == b"ATTTACTATGCCATGCTGCTGCT"


@pytest.mark.django_db
def test_primer_28727_right_head(primer_28727_right):
    assert primer_28727_right.head.seq == b"GGCATCTGTTCTTGTTTGCACAT"


@pytest.mark.django_db
def test_primer_28727_left_tail(primer_28727_left):
    assert primer_28727_left.tail.seq == b"CTACACGACGCTCTTCCGATCT"


@pytest.mark.django_db
def test_primer_28727_right_tail(primer_28727_right):
    assert primer_28727_right.tail.seq == b"CAGACGTGTGCTCTTCCGATCT"


@pytest.mark.django_db
def test_primer_28727_left(primer_28727_left):
    assert primer_28727_left.sequence.seq == b"CTACACGACGCTCTTCCGATCTATTTACTATGCCATGCTGCTGCT"


@pytest.mark.django_db
def test_primer_28727_right(primer_28727_right):
    assert primer_28727_right.sequence.seq == b"CAGACGTGTGCTCTTCCGATCTGGCATCTGTTCTTGTTTGCACAT"


@pytest.mark.django_db
def test_primer_28734_left(primer_28734_left):
    assert primer_28734_left.sequence.seq == b"CTACACGACGCTCTTCCGATCTAAAGGCTTCTCCCCATTCCAAAG"


@pytest.mark.django_db
def test_primer_28734_right(primer_28734_right):
    assert primer_28734_right.sequence.seq == b"CAGACGTGTGCTCTTCCGATCTGGAAGTAGTGTGTGCTTGGACTT"
