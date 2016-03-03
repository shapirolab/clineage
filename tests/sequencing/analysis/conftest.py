import pytest
import os
import uuid
from django.conf import settings

from sequencing.analysis.models import SampleReads, AdamMergedReads, \
    AdamReadsIndex, AdamMarginAssignment, AdamAmpliconReads
from sequencing.analysis.models import LEFT, RIGHT
from misc.utils import get_unique_path

from tests.sequencing.runs.conftest import *
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.unwrapping.conftest import *

file_fixtures_path = os.path.join(*(os.path.split(os.path.dirname(os.path.realpath(__file__)))[:-1] + ("ngs_fixtures",)))


@pytest.fixture()
def samplereads(demultiplexing, magicalpcr1barcodedcontent, magicalpcr1library):
    fastq_r1 = get_unique_path("fastq")
    fastq_r2 = get_unique_path("fastq")
    os.symlink(
        os.path.join(file_fixtures_path, '1448-Viktor-AAR20-BC81_S321_L001_R1_001/28727_and_28734_R1.fastq'),
        fastq_r1
    )
    os.symlink(
        os.path.join(file_fixtures_path, '1448-Viktor-AAR20-BC81_S321_L001_R1_001/28727_and_28734_R2.fastq'),
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
def merged_reads_stripped_fasta():
    return {
        (
            'M00393:167:000000000-ABF3N:1:1101:20583:16769',
            'AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1101:7043:8470',
            'AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1102:6086:12380',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1106:14231:8476',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1109:4497:5075',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1112:13016:19696',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1112:20425:16124',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2101:12593:17954',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2104:19595:12515',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2109:21650:4662',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2115:26810:9430',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2117:19711:1998',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ),
    }


@pytest.fixture()
def mergedreads(samplereads):
    src_prefix = os.path.join(file_fixtures_path,
                              '1448-Viktor-AAR20-BC81_S321_L001_R1_001/1')
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
        demux_reads=samplereads,
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
                           '1448-Viktor-AAR20-BC81_S321_L001_R1_001/ind1M2')
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
                           '1448-Viktor-AAR20-BC81_S321_L001_R1_001/ind1F')
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
def adammarginassignment(readsindex_fwd_and_merged, require_unwrappers):
    alignment_file_name = get_unique_path("sam")
    os.symlink(
        os.path.join(file_fixtures_path,
                     '1448-Viktor-AAR20-BC81_S321_L001_R1_001/margine_assignemnt.sam'),
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
def reads_unwrappers(pu_28734):
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
def reads_by_unwrappers(pu_28734):
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
                     '1448-Viktor-AAR20-BC81_S321_L001_R1_001/28734.fastq'),
        fastq_path
    )
    aar = AdamAmpliconReads.objects.create(
        margin_assignment=adammarginassignment,
        unwrapper=pu_28734,
        fastq=fastq_path
    )
    return aar
