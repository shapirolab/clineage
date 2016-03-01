import pytest

from sequencing.analysis.models import DemultiplexedReads, MergingScheme, MergedReads, ReadsIndex, \
    UGSAssignment

from tests.sequencing.runs.conftest import *
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.unwrapping.conftest import *


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
    prefix = '/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/1'
    mr = MergedReads.objects.create(
        demux_reads=demultiplexedreads,
        merging_scheme=mergingscheme,
        assembled_fastq="{}.assembled.fastq".format(prefix),
        discarded_fastq="{}.discarded.fastq".format(prefix),
        unassembled_forward_fastq="{}.unassembled.forward.fastq".format(prefix),
        unassembled_reverse_fastq="{}.unassembled.reverse.fastq".format(prefix),
    )
    return mr


@pytest.fixture()
def readsindex_merged_only(mergedreads):
    ri = ReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='M',  # Only merged
        index_dump_dir='/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/ind1M2',
        padding=5,
    )
    return ri


@pytest.fixture()
def readsindex_fwd_and_merged(mergedreads):
    ri = ReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='F',  # Merged and unassembled_forward
        index_dump_dir='/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/ind1F',
        padding=5,
    )
    return ri


@pytest.fixture()
def ugsassignment(readsindex_fwd_and_merged):
    ugsa = UGSAssignment.objects.create(
        reads_index=readsindex_fwd_and_merged,
    )
    return ugsa




