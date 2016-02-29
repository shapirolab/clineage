import pytest

from sequencing.analysis.models import DemultiplexingScheme, Demultiplexing, DemultiplexedReads, MergingScheme, MergedReads, ReadsIndex, \
    UGSAssignment

from tests.sequencing.runs.conftest import *
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.unwrapping.conftest import *


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


@pytest.fixture()
def readsindex_merged_only(mergedreads):
    ri = ReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='M',  # Merged and unassembled_forward
        padding=5,
    )
    return ri


@pytest.fixture()
def readsindex_fwd_and_merged(mergedreads):
    ri = ReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='F',  # Merged and unassembled_forward
        padding=5,
    )
    return ri


@pytest.fixture()
def ugsassignment(readsindex_fwd_and_merged):
    ugsa = UGSAssignment.objects.create(
        reads_index=readsindex_fwd_and_merged,
    )
    return ugsa




