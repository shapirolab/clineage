import pytest


@pytest.mark.django_db
def test_ter(ter_28727):
    assert ter_28727.left_primer.sequence.seq == b"CTACACGACGCTCTTCCGATCTATTTACTATGCCATGCTGCTGCT"
    assert "{}".format(ter_28727) == "chrX_81316201_81316236_ACG_fwd, Seq19272_rev"
