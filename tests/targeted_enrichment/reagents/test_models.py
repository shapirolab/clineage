import pytest


@pytest.mark.django_db
def test_ter(ter_28727):
    assert ter_28727.left_primer.sequence.seq == b"CTACACGACGCTCTTCCGATCTATTTACTATGCCATGCTGCTGCT"
    assert ter_28727.te.targets.get(name='X_81316201_81316236').slice.start_pos == 81316201
