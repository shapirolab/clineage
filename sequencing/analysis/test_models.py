import pytest
import os
from models import DemultiplexingScheme, Demultiplexing, DemultiplexedReads, MergingScheme, MergedReads

from accounts.test_user import user
from misc.test_models import human_taxa
from linapp.test_models import protocoltype
from primers.parts.test_models import illuminareadingadaptor1, illuminareadingadaptor2, \
    illuminareadingadaptor1cuts, illuminareadingadaptor2cuts, \
    dnabarcode1, dnabarcode2, dnabarcode1_a, dnabarcode2_a
from sampling.test_models import human_cell, human_individual, composition
from lib_prep.workflows.test_models import barcodepair, barcodepair_a, amplifiedcontent, cellcontentprotocol
from lib_prep.multiplexes.test_models import panel, pcr1multiplex, pcr1multiplexcollection
from ..runs.test_models import machine, ngskit, machinetype
from genomes.test_models import hg19_assembly, hg19_chromosome, \
    slice_28727_left, slice_28727_right, slice_28727_target_a, slice_28727_target_b,\
    slice_28734_left, slice_28734_right, slice_28734_target_a
from targeted_enrichment.planning.test_models import ugs_28727_left, ugs_28727_right, \
    ugs_28734_left, ugs_28734_right, ms_28727_a, ms_28727_b, ms_28734_a
from primers.synthesis.test_models import primer_28727_left, primer_28727_right, \
    primer_28734_left, primer_28734_right
from targeted_enrichment.planning.test_models import te_28727, te_28734
from targeted_enrichment.reagents.test_models import ter_28727, ter_28734

from ..runs.test_models import ngsrun
from lib_prep.workflows.test_models import magicalpcr1library, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a


@pytest.fixture()
def demultiplexingscheme(db):
    ds = DemultiplexingScheme.objects.create(
        name='test demux scheme',
        description='wrovhnwpovnwecpqkewmc',
    )
    return ds


@pytest.fixture()
def demultiplexing(demultiplexingscheme, ngsrun):
    dm = Demultiplexing.objects.create(
        ngs_run=ngsrun,
        demux_scheme=demultiplexingscheme
    )
    return dm


@pytest.fixture()
def demultiplexedreads(demultiplexing, magicalpcr1barcodedcontent, magicalpcr1library):
    dr = DemultiplexedReads.objects.create(
        demux=demultiplexing,
        barcoded_content=magicalpcr1barcodedcontent,
        library=magicalpcr1library,
        fastq1='/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/28727_and_28734_R1.fastq',
        fastq2='/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/28727_and_28734_R2.fastq',
    )
    return dr


@pytest.fixture()
def mergingscheme(db):
    ms = MergingScheme.objects.create(
        name='test merging scheme',
        description='sdnvjweivobwvciwenc wsnfcueqwlcnewqc',
    )
    return ms


@pytest.fixture()
def mergedreads(demultiplexedreads, mergingscheme):
    mr = MergedReads.objects.create(
        demux_read=demultiplexedreads,
        merge_scheme=mergingscheme,
    )
    return mr


@pytest.mark.django_db
def test_demultiplexedreads(demultiplexedreads):
    assert demultiplexedreads.barcoded_content.content.name == 'human amplified content'


@pytest.mark.django_db
def test_mergedreads(mergedreads):
    assert os.path.isfile(mergedreads.demux_read.fastq1)
    assert os.path.isfile(mergedreads.demux_read.fastq2)


@pytest.mark.django_db
def test_runmerge(mergedreads):
    mergedreads.run_merge()
    assert os.path.isfile(mergedreads.assembled_fastq)
    assert os.path.isfile(mergedreads.discarded_fastq)
    assert os.path.isfile(mergedreads.unassembled_forward_fastq)
    assert os.path.isfile(mergedreads.unassembled_reverse_fastq)