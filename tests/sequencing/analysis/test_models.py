import pytest
import os
from Bio import SeqIO


def strip_fasta_records(fasta_records):
    for rec in fasta_records:
        yield rec.id, str(rec.seq)


@pytest.mark.django_db
def test_samplereads(samplereads):
    assert samplereads.barcoded_content.content.name == 'human amplified content'


@pytest.mark.django_db
def test_mergedreads(mergedreads):
    assert os.path.isfile(mergedreads.demux_reads.fastq1)
    assert os.path.isfile(mergedreads.demux_reads.fastq2)
    assert set(strip_fasta_records(mergedreads.included_reads_generator('M'))) == {
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
    assert set(strip_fasta_records(mergedreads.included_reads_generator('M'))) == set(strip_fasta_records(mergedreads.included_reads_generator('F')))


@pytest.mark.django_db
def test_runmerge(samplereads, mergingscheme):
    mr = samplereads.run_merge(mergingscheme)
    assert mr.demux_reads is samplereads
    assert os.path.isfile(mr.assembled_fastq)
    assert os.path.isfile(mr.discarded_fastq)
    assert os.path.isfile(mr.unassembled_forward_fastq)
    assert os.path.isfile(mr.unassembled_reverse_fastq)
    mr.delete()
    assert not os.path.exists(mr.assembled_fastq)
    assert not os.path.exists(mr.discarded_fastq)
    assert not os.path.exists(mr.unassembled_forward_fastq)
    assert not os.path.exists(mr.unassembled_reverse_fastq)


# @pytest.mark.django_db
# def test_readsindex_merged_only(readsindex_merged_only):
    # assert os.path.isfile(readsindex_merged_only.merged_reads.assembled_fastq)
    # readsindex_merged_only.create_final_merged_fastq()
    # assert os.path.isfile(readsindex_merged_only.padded_reads_fasta)
# 
# 
# @pytest.mark.django_db
# def test_readsindex_fwd_and_merged(readsindex_fwd_and_merged):
    # assert os.path.isfile(readsindex_fwd_and_merged.merged_reads.assembled_fastq)
    # assert os.path.isfile(readsindex_fwd_and_merged.merged_reads.unassembled_forward_fastq)
    # readsindex_fwd_and_merged.create_final_merged_fastq()
    # assert os.path.isfile(readsindex_fwd_and_merged.padded_reads_fasta)


@pytest.mark.django_db
def test_readsindex_bowtie2build_merged_only(mergedreads):
    assert os.path.isfile(mergedreads.assembled_fastq)
    ri = mergedreads.create_reads_index(included_reads='M',padding=5)
    assert set(os.path.join(ri.index_dump_dir, x) for x in os.listdir(ri.index_dump_dir)) == set(ri.files)
    ri.delete()
    assert not os.path.exists(ri.index_dump_dir)


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
