import pytest


@pytest.mark.django_db
def test_pcr1multiplex(pcr1multiplex, te_28727):
    assert pcr1multiplex.name == 'test PCR1Multiplex'
    assert pcr1multiplex.ters.count() == 8
    assert pcr1multiplex.ters.select_subclasses().get(te=te_28727)\
        .left_primer.sequence.seq == b"CTACACGACGCTCTTCCGATCTATTTACTATGCCATGCTGCTGCT"


@pytest.mark.django_db
def test_pcr1panel(pcr1panel, te_28727):
    assert pcr1panel.name == 'test PCR1Panel'
    assert pcr1panel.mpxs.count() == 1
    assert pcr1panel.mpxs.get().ters.count() == 8
