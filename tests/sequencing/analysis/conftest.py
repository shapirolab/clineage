import pytest
import os
import uuid
from django.conf import settings

from sequencing.analysis.models import DemultiplexedReads, MergingScheme, MergedReads, ReadsIndex, \
    UGSAssignment

from misc.utils import get_unique_path

from tests.sequencing.runs.conftest import *
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.unwrapping.conftest import *


@pytest.fixture()
def demultiplexedreads(demultiplexing, magicalpcr1barcodedcontent, magicalpcr1library):
    fastq_r1 = get_unique_path("fastq")
    fastq_r2 = get_unique_path("fastq")
    os.symlink(
        '/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/28727_and_28734_R1.fastq',
        fastq_r1
    )
    os.symlink(
        '/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/28727_and_28734_R2.fastq',
        fastq_r2
    )
    dr = DemultiplexedReads.objects.create(
        demux=demultiplexing,
        barcoded_content=magicalpcr1barcodedcontent,
        library=magicalpcr1library,
        fastq1=fastq_r1,
        fastq2=fastq_r2,
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
    src_prefix = '/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/1'
    dst_prefix = get_unique_path()
    os.symlink(
        "{}.assembled.fastq".format(src_prefix),
        "{}.assembled.fastq".format(dst_prefix)
    )
    os.symlink(
        "{}.discarded.fastq".format(src_prefix),
        "{}.discarded.fastq".format(dst_prefix)
    )
    os.symlink(
        "{}.unassembled.forward.fastq".format(src_prefix),
        "{}.unassembled.forward.fastq".format(dst_prefix)
    )
    os.symlink(
        "{}.unassembled.reverse.fastq".format(src_prefix),
        "{}.unassembled.reverse.fastq".format(dst_prefix)
    )
    mr = MergedReads.objects.create(
        demux_reads=demultiplexedreads,
        merging_scheme=mergingscheme,
        assembled_fastq="{}.assembled.fastq".format(dst_prefix),
        discarded_fastq="{}.discarded.fastq".format(dst_prefix),
        unassembled_forward_fastq="{}.unassembled.forward.fastq".format(dst_prefix),
        unassembled_reverse_fastq="{}.unassembled.reverse.fastq".format(dst_prefix),
    )
    return mr


def link_index_files(src_dir, dst_dir):
    os.mkdir(dst_dir)
    for index_file in os.listdir(src_dir):
        os.symlink(
            os.path.join(src_dir, index_file),
            os.path.join(dst_dir, index_file),
        )


@pytest.fixture()
def readsindex_merged_only(mergedreads):
    src_dir = '/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/ind1M2'
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    ri = ReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='M',  # Only merged
        index_dump_dir=dst_dir,
        padding=5,
    )
    return ri


@pytest.fixture()
def readsindex_fwd_and_merged(mergedreads):
    src_dir = '/net/mraid11/export/dcstor/Ofir/ngs_fixtures/1448-Viktor-AAR20-BC81_S321_L001_R1_001/ind1F'
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    ri = ReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='F',  # Merged and unassembled_forward
        index_dump_dir=dst_dir,
        padding=5,
    )
    return ri


@pytest.fixture()
def ugsassignment(readsindex_fwd_and_merged):
    ugsa = UGSAssignment.objects.create(
        reads_index=readsindex_fwd_and_merged,
    )
    return ugsa




