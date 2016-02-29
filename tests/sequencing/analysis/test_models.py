import pytest
import os
from Bio import SeqIO


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


@pytest.mark.django_db
def test_readsindex_merged_only(readsindex_merged_only):
    readsindex_merged_only.merged_reads.run_merge()
    assert os.path.isfile(readsindex_merged_only.merged_reads.assembled_fastq)
    readsindex_merged_only.create_final_merged_fastq()
    assert os.path.isfile(readsindex_merged_only.padded_reads_fasta)


@pytest.mark.django_db
def test_readsindex_fwd_and_merged(readsindex_fwd_and_merged):
    readsindex_fwd_and_merged.merged_reads.run_merge()
    assert os.path.isfile(readsindex_fwd_and_merged.merged_reads.assembled_fastq)
    assert os.path.isfile(readsindex_fwd_and_merged.merged_reads.unassembled_forward_fastq)
    readsindex_fwd_and_merged.create_final_merged_fastq()
    assert os.path.isfile(readsindex_fwd_and_merged.padded_reads_fasta)


@pytest.mark.django_db
def test_readsindex_bowtie2build(readsindex_merged_only):
    readsindex_merged_only.merged_reads.run_merge()
    assert os.path.isfile(readsindex_merged_only.merged_reads.assembled_fastq)
    readsindex_merged_only.create_final_merged_fastq()
    assert os.path.isfile(readsindex_merged_only.padded_reads_fasta)
    readsindex_merged_only.index_reads()
    expected_files = []
    for i in xrange(1, 5):
        f = "{}.{}.bt2".format(readsindex_merged_only.index_files_prefix, i)
        assert os.path.isfile(f)
        expected_files.append(f)
    for i in xrange(1, 3):
        f = "{}.rev.{}.bt2".format(readsindex_merged_only.index_files_prefix, i)
        assert os.path.isfile(f)
        expected_files.append(f)
    for f in readsindex_merged_only.collect_bt_files():
        if f not in expected_files:
            assert False


def strip_fasta_records(fasta_records):
    for rec in fasta_records:
        yield rec.id, str(rec.seq)


@pytest.mark.django_db
def test_ugsassignment_create_panel_fasta(ugsassignment, pu_28727, pu_28734):
    ugsassignment.create_panel_fasta()
    assert os.path.isfile(ugsassignment.panel_fasta)
    assert set(strip_fasta_records(SeqIO.parse(ugsassignment.panel_fasta, "fasta"))) == {
        ("unwrapper_1_left", "TTTACTATGCCATGCTGCTGCT",),
        ("unwrapper_1_right", "TGTGCAAACAAGAACAGATGCC",),
        ("unwrapper_2_left", "AAGGCTTCTCCCCATTCCAAAG",),
        ("unwrapper_2_right", "AGTCCAAGCACACACTACTTCC",),
    }


@pytest.mark.django_db
def test_ugsassignment_align_primers_to_reads_basic(ugsassignment, pu_28727, pu_28734):
    ugsassignment.create_panel_fasta()
    assert os.path.isfile(ugsassignment.panel_fasta)
    ugsassignment.align_primers_to_reads()
    assert os.path.isfile(ugsassignment.primer_reads_alignment)


@pytest.mark.django_db
def test_ugsassignment_align_primers_to_reads_collection(ugsassignment, pu_28727, pu_28734):
    ugsassignment.create_panel_fasta()
    assert os.path.isfile(ugsassignment.panel_fasta)
    ugsassignment.align_primers_to_reads()
    assert os.path.isfile(ugsassignment.primer_reads_alignment)
    assert ugsassignment.collect_mappings_from_sam() == {
        'unwrapper_2_right': {
            'M00393:167:000000000-ABF3N:1:2115:26810:9430',
            'M00393:167:000000000-ABF3N:1:2101:12593:17954',
            'M00393:167:000000000-ABF3N:1:1112:13016:19696',
            'M00393:167:000000000-ABF3N:1:1101:7043:8470',
            'M00393:167:000000000-ABF3N:1:1112:20425:16124',
            'M00393:167:000000000-ABF3N:1:1101:20583:16769',
            'M00393:167:000000000-ABF3N:1:2117:19711:1998',
            'M00393:167:000000000-ABF3N:1:1102:6086:12380',
            'M00393:167:000000000-ABF3N:1:1109:4497:5075',
            'M00393:167:000000000-ABF3N:1:2109:21650:4662',
            'M00393:167:000000000-ABF3N:1:1106:14231:8476',
            'M00393:167:000000000-ABF3N:1:2104:19595:12515'
        },
        'unwrapper_2_left': {
            'M00393:167:000000000-ABF3N:1:2117:19711:1998',
            'M00393:167:000000000-ABF3N:1:1106:14231:8476',
            'M00393:167:000000000-ABF3N:1:2115:26810:9430',
            'M00393:167:000000000-ABF3N:1:2104:19595:12515',
            'M00393:167:000000000-ABF3N:1:2109:21650:4662',
            'M00393:167:000000000-ABF3N:1:1112:13016:19696',
            'M00393:167:000000000-ABF3N:1:1109:4497:5075',
            'M00393:167:000000000-ABF3N:1:2101:12593:17954',
            'M00393:167:000000000-ABF3N:1:1102:6086:12380',
            'M00393:167:000000000-ABF3N:1:1112:20425:16124',
            'M00393:167:000000000-ABF3N:1:1101:7043:8470',
            'M00393:167:000000000-ABF3N:1:1101:20583:16769'
        }
    }


@pytest.mark.django_db
def test_ugsassignment_read_ids_by_te(ugsassignment, pu_28727, pu_28734):
    ugsassignment.create_panel_fasta()
    ugsassignment.align_primers_to_reads()
    assert os.path.isfile(ugsassignment.primer_reads_alignment)
    assert list(ugsassignment.read_ids_by_unwrapper()) == [
        (
            pu_28727,
            set()
        ),
        (
            pu_28734,
            {
                'M00393:167:000000000-ABF3N:1:2115:26810:9430',
                'M00393:167:000000000-ABF3N:1:2101:12593:17954',
                'M00393:167:000000000-ABF3N:1:1112:13016:19696',
                'M00393:167:000000000-ABF3N:1:1101:7043:8470',
                'M00393:167:000000000-ABF3N:1:1112:20425:16124',
                'M00393:167:000000000-ABF3N:1:1101:20583:16769',
                'M00393:167:000000000-ABF3N:1:2117:19711:1998',
                'M00393:167:000000000-ABF3N:1:1102:6086:12380',
                'M00393:167:000000000-ABF3N:1:1109:4497:5075',
                'M00393:167:000000000-ABF3N:1:2109:21650:4662',
                'M00393:167:000000000-ABF3N:1:1106:14231:8476',
                'M00393:167:000000000-ABF3N:1:2104:19595:12515'
            }
        )
    ]


@pytest.mark.django_db
def test_ugsassignment_reads_by_te(ugsassignment, pu_28727, pu_28734):
    ugsassignment.reads_index.merged_reads.run_merge()
    ugsassignment.create_panel_fasta()
    ugsassignment.align_primers_to_reads()
    assert os.path.isfile(ugsassignment.primer_reads_alignment)
    uws_and_reads = list(ugsassignment.reads_by_unwrapper())
    uw, reads_list = uws_and_reads[0]
    assert uw == pu_28727
    assert set(strip_fasta_records(reads_list)) == set()
    uw, reads_list = uws_and_reads[1]
    assert uw == pu_28734
    assert set(strip_fasta_records(reads_list)) == {
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
