import pytest


@pytest.mark.django_db
def test_pcr1multiplex(pcr1multiplex, te_28727):
    assert pcr1multiplex.name == 'test PCR1Multiplex'
    assert pcr1multiplex.ters.count() == 2
    assert pcr1multiplex.ters.select_subclasses().get(te=te_28727)\
        .left_primer.sequence.seq == "CTACACGACGCTCTTCCGATCTATTTACTATGCCATGCTGCTGCT"
    assert pcr1multiplex.ters.select_subclasses().get(te=te_28727)\
        .te.targets.get(name='X_81316201_81316236').slice.start_pos == 81316201


@pytest.mark.django_db
def test_pcr1multiplexcollection(pcr1multiplexcollection, te_28727):
    assert pcr1multiplexcollection.panel.name == 'test Panel'
    assert pcr1multiplexcollection.mpxs.count() == 1
    assert pcr1multiplexcollection.mpxs.get().ters.count() == 2
    assert pcr1multiplexcollection.mpxs.get().ters.select_subclasses().get(te=te_28727)\
        .te.targets.get(name='X_81316201_81316236').slice.start_pos == 81316201
