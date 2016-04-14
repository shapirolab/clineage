import pytest


@pytest.mark.django_db
def test_pu_28727(pu_28727):
    assert pu_28727.left_margin.seq == b"ATTTACTATGCCATGCTGCTGCT"
