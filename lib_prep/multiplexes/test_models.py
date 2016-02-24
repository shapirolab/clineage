import pytest
from models import PCR1Multiplex, Panel, PCR1MultiplexCollection

from misc.test_models import human_taxa
from genomes.test_models import hg19_assembly, hg19_chromosome, \
    slice_28727_left, slice_28727_right, slice_28727_target_a, slice_28727_target_b,\
    slice_28734_left, slice_28734_right, slice_28734_target_a
from targeted_enrichment.planning.test_models import ugs_28727_left, ugs_28727_right, \
    ugs_28734_left, ugs_28734_right, ms_28727_a, ms_28727_b, ms_28734_a
from primers.synthesis.test_models import primer_28727_left, primer_28727_right, \
    primer_28734_left, primer_28734_right
from primers.parts.test_models import illuminareadingadaptor1cuts, illuminareadingadaptor2cuts, \
    illuminareadingadaptor1, illuminareadingadaptor2

from targeted_enrichment.planning.test_models import te_28727, te_28734
from targeted_enrichment.reagents.test_models import ter_28727, ter_28734


@pytest.fixture()
def pcr1multiplex(ter_28727, ter_28734):
    pcr1m = PCR1Multiplex.objects.create(
        name='test PCR1Multiplex'
    )
    pcr1m.ters = [ter_28727, ter_28734]
    return pcr1m


@pytest.fixture()
def panel(te_28727, te_28734):
    p = Panel.objects.create(
        name='test Panel'
    )
    p.tes = [te_28727, te_28734]
    return p


@pytest.fixture()
def pcr1multiplexcollection(panel, pcr1multiplex):
    pcr1mc = PCR1MultiplexCollection.objects.create(
        panel=panel
    )
    pcr1mc.mpxs = [pcr1multiplex]
    return pcr1mc


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
