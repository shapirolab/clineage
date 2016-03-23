import pytest
import os
import uuid
from django.conf import settings
from Bio import SeqIO

from sequencing.analysis.models import SampleReads, AdamMergedReads, \
    AdamReadsIndex, AdamMarginAssignment, AdamAmpliconReads, \
    AdamMSVariations, AdamHistogram
from sequencing.analysis.models import LEFT, RIGHT
from misc.utils import get_unique_path

from tests.sequencing.runs.conftest import *
from tests.sequencing.analysis.pu_28727_adam_ms_variations import VARS_28727
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.amplicons.conftest import *
from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM
from tests.flat_dict import FlatDict
from tests.sequencing.analysis.reads_dict import READS_DICT_ADAM, ASSEMBLED, \
    UNASSEMBLED


file_fixtures_path = os.path.join(*(os.path.split(os.path.dirname(os.path.realpath(__file__)))[:-1] + ("ngs_fixtures",)))


@pytest.fixture(scope="session")
def adam_reads_fd():
    return FlatDict(READS_DICT_ADAM, [R1, R2, RM])


@pytest.yield_fixture(scope="session")
def sample_reads_files_d(adam_reads_fd):
    d = {}
    for k, r_d in adam_reads_fd.items():
        fastq_r1 = get_unique_path("fastq")
        fastq_r2 = get_unique_path("fastq")
        SeqIO.write(r_d[R1], fastq_r1, "fastq")
        SeqIO.write(r_d[R2], fastq_r2, "fastq")
        d[k] = {R1: fastq_r1, R2: fastq_r2}
    yield d
    for f_d in d.itervalues():
        os.unlink(f_d[R1])
        os.unlink(f_d[R2])


@pytest.yield_fixture()
def sample_reads_d(sample_reads_files_d, demultiplexing, magicalpcr1barcodedcontent, magicalpcr1library):
    d = {}
    for k, f_d in sample_reads_files_d.iteritems():
        fastq_r1 = get_unique_path("fastq")
        fastq_r2 = get_unique_path("fastq")
        os.symlink(f_d[R1], fastq_r1)
        os.symlink(f_d[R2], fastq_r2)
        sr = SampleReads.objects.create(
            demux=demultiplexing,
            barcoded_content=magicalpcr1barcodedcontent,
            library=magicalpcr1library,
            fastq1=fastq_r1,
            fastq2=fastq_r2,
        )
        d[k] = sr
    yield d
    for sr in d.itervalues():
        # sr.delete()
        os.unlink(sr.fastq1)
        os.unlink(sr.fastq2)


@pytest.yield_fixture(scope="session")
def adam_merged_reads_files_d(adam_reads_fd):
    d = {}
    for k in adam_reads_fd.keys():
        assembled_fastq = get_unique_path("fastq")
        unassembled_forward_fastq = get_unique_path("fastq")
        unassembled_reverse_fastq = get_unique_path("fastq")
        discarded_fastq = get_unique_path("fastq")
        SeqIO.write(adam_reads_fd[k, ASSEMBLED][RM], assembled_fastq, "fastq")
        SeqIO.write(adam_reads_fd[k, UNASSEMBLED][R1], unassembled_forward_fastq, "fastq")
        SeqIO.write(adam_reads_fd[k, UNASSEMBLED][R2], unassembled_reverse_fastq, "fastq")
        SeqIO.write((), discarded_fastq, "fastq")
        d[k] = {
            ASSEMBLED: assembled_fastq,
            (UNASSEMBLED, R1): unassembled_forward_fastq,
            (UNASSEMBLED, R2): unassembled_reverse_fastq,
            None: discarded_fastq
        }
    yield d
    for f_d in d.itervalues():
        os.unlink(f_d[ASSEMBLED])
        os.unlink(f_d[UNASSEMBLED, R1])
        os.unlink(f_d[UNASSEMBLED, R2])
        os.unlink(f_d[None])


@pytest.yield_fixture()
def adam_merged_reads_d(adam_merged_reads_files_d, sample_reads_d):
    d = {}
    for k, f_d in adam_merged_reads_files_d.iteritems():
        dst_prefix = get_unique_path()
        os.symlink(
            f_d[ASSEMBLED],
            "{}.assembled.fastq".format(dst_prefix)
        )
        os.symlink(
            f_d[None],
            "{}.discarded.fastq".format(dst_prefix)
        )
        os.symlink(
            f_d[UNASSEMBLED, R1],
            "{}.unassembled.forward.fastq".format(dst_prefix)
        )
        os.symlink(
            f_d[UNASSEMBLED, R2],
            "{}.unassembled.reverse.fastq".format(dst_prefix)
        )
        mr = AdamMergedReads.objects.create(
            sample_reads=sample_reads_d[k],
            assembled_fastq="{}.assembled.fastq".format(dst_prefix),
            discarded_fastq="{}.discarded.fastq".format(dst_prefix),
            unassembled_forward_fastq="{}.unassembled.forward.fastq".format(dst_prefix),
            unassembled_reverse_fastq="{}.unassembled.reverse.fastq".format(dst_prefix),
        )
        d[k] = mr
    yield d
    for mr in d.itervalues():
        # mr.delete()
        os.unlink(mr.assembled_fastq)
        os.unlink(mr.discarded_fastq)
        os.unlink(mr.unassembled_forward_fastq)
        os.unlink(mr.unassembled_reverse_fastq)


@pytest.fixture()
def samplereads(demultiplexing, magicalpcr1barcodedcontent, magicalpcr1library):
    fastq_r1 = get_unique_path("fastq")
    fastq_r2 = get_unique_path("fastq")
    os.symlink(
        os.path.join(file_fixtures_path, 'adam/bcs/bc1/bc1_R1.fastq'),
        fastq_r1
    )
    os.symlink(
        os.path.join(file_fixtures_path, 'adam/bcs/bc1/bc1_R2.fastq'),
        fastq_r2
    )
    sr = SampleReads.objects.create(
        demux=demultiplexing,
        barcoded_content=magicalpcr1barcodedcontent,
        library=magicalpcr1library,
        fastq1=fastq_r1,
        fastq2=fastq_r2,
    )
    return sr


@pytest.fixture()
def mergedreads(samplereads):
    src_prefix = os.path.join(file_fixtures_path,
                              'adam/bcs/bc1/bc1')
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
    mr = AdamMergedReads.objects.create(
        sample_reads=samplereads,
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
    src_dir = os.path.join(file_fixtures_path,
                           'adam/bcs/bc1/M/')
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    ri = AdamReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='M',  # Only merged
        index_dump_dir=dst_dir,
        padding=5,
    )
    return ri


@pytest.fixture()
def readsindex_fwd_and_merged(mergedreads):
    src_dir = os.path.join(file_fixtures_path,
                           'adam/bcs/bc1/F/')
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    ri = AdamReadsIndex.objects.create(
        merged_reads=mergedreads,
        included_reads='F',  # Merged and unassembled_forward
        index_dump_dir=dst_dir,
        padding=5,
    )
    return ri


@pytest.fixture()
def adammarginassignment(readsindex_fwd_and_merged, require_amplicons):
    alignment_file_name = get_unique_path("sam")
    os.symlink(
        os.path.join(file_fixtures_path,
                     'adam/bcs/bc1/F/margins.sam'),
        alignment_file_name
    )
    ama = AdamMarginAssignment.objects.create(
        reads_index=readsindex_fwd_and_merged,
        assignment_sam=alignment_file_name,
    )
    return ama


@pytest.fixture()
def reads_matches(pu_28734):
    return {
        'M00393:167:000000000-ABF3N:1:1101:20583:16769':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1101:7043:8470':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1102:6086:12380':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:1106:14231:8476':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:1109:4497:5075':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1112:13016:19696':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1112:20425:16124':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:2101:12593:17954':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:2104:19595:12515':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:2109:21650:4662':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:2115:26810:9430':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:2117:19711:1998':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)}
    }


@pytest.fixture()
def reads_amplicons(pu_28734):
    return {
        ('M00393:167:000000000-ABF3N:1:1101:20583:16769',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1101:7043:8470',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1102:6086:12380',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1106:14231:8476',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1109:4497:5075',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1112:13016:19696',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1112:20425:16124',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2101:12593:17954',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2104:19595:12515',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2109:21650:4662',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2115:26810:9430',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2117:19711:1998',
            pu_28734),
    }


@pytest.fixture()
def reads_by_amplicons(pu_28734):
    return {
        pu_28734: {
            'M00393:167:000000000-ABF3N:1:1101:20583:16769',
            'M00393:167:000000000-ABF3N:1:1101:7043:8470',
            'M00393:167:000000000-ABF3N:1:1102:6086:12380',
            'M00393:167:000000000-ABF3N:1:1106:14231:8476',
            'M00393:167:000000000-ABF3N:1:1109:4497:5075',
            'M00393:167:000000000-ABF3N:1:1112:13016:19696',
            'M00393:167:000000000-ABF3N:1:1112:20425:16124',
            'M00393:167:000000000-ABF3N:1:2101:12593:17954',
            'M00393:167:000000000-ABF3N:1:2104:19595:12515',
            'M00393:167:000000000-ABF3N:1:2109:21650:4662',
            'M00393:167:000000000-ABF3N:1:2115:26810:9430',
            'M00393:167:000000000-ABF3N:1:2117:19711:1998',
        }
    }


@pytest.fixture()
def adamampliconreads(adammarginassignment, pu_28734):
    fastq_path = get_unique_path("fastq")
    os.symlink(
        os.path.join(file_fixtures_path,
                     'adam/bcs/bc1/F/28734/28734.fastq'),
        fastq_path
    )
    fastq1_path = get_unique_path("fastq")
    os.symlink(
        os.path.join(file_fixtures_path,
                     'adam/bcs/bc1/F/28734/28734_R1.fastq'),
        fastq1_path
    )
    fastq2_path = get_unique_path("fastq")
    os.symlink(
        os.path.join(file_fixtures_path,
                     'adam/bcs/bc1/F/28734/28734_R2.fastq'),
        fastq2_path
    )
    aar = AdamAmpliconReads.objects.create(
        margin_assignment=adammarginassignment,
        amplicon=pu_28734,
        fastq=fastq_path,
        fastq1=fastq1_path,
        fastq2=fastq2_path,
    )
    return aar


@pytest.fixture()
def pu_28727_adam_ms_variations():
    return VARS_28727


@pytest.fixture()
def amsv_28727(pu_28727):
    src_dir = os.path.join(file_fixtures_path, 'amsv_28727')
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    amsv = AdamMSVariations.objects.create(
        amplicon=pu_28727,
        index_dump_dir=dst_dir,
        padding=50,
    )
    return amsv
    

@pytest.fixture()
def amsv_28734(pu_28734):
    src_dir = os.path.join(file_fixtures_path, 'amsv_28734')
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    amsv = AdamMSVariations.objects.create(
        amplicon=pu_28734,
        index_dump_dir=dst_dir,
        padding=50,
    )
    return amsv
    

@pytest.fixture()
def require_adammsvariations(amsv_28727, amsv_28734):
    pass


@pytest.fixture()
def adamhistogram(adamampliconreads):
    alignment_file_name = get_unique_path("sam")
    os.symlink(
        os.path.join(file_fixtures_path,
                     'adam/bcs/bc1/F/28734/28734.sam'),
        alignment_file_name
    )
    ama = AdamHistogram.objects.create(
        sample_reads=adamampliconreads.margin_assignment.reads_index \
            .merged_reads.sample_reads,
        amplicon_reads=adamampliconreads,
        assignment_sam=alignment_file_name,
    )
    return ama
