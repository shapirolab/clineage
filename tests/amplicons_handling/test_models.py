import pytest


@pytest.mark.django_db
def ugs_28734(primer_28734_left):
    assert primer_28734_left.sequence.seq == b"CTACACGACGCTCTTCCGATCTAAAGGCTTCTCCCCATTCCAAAG"