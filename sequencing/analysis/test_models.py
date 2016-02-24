import pytest
import os
from Bio import SeqIO
from models import DemultiplexingScheme, Demultiplexing, DemultiplexedReads, MergingScheme, MergedReads, ReadsIndex, \
    UGSAssignment

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
from targeted_enrichment.planning.test_models import ugs_28727_left, ugs_28727_right, ugs_28734_left, ugs_28734_right, \
    te_28727, te_28734, ms_28727_a, ms_28727_b, ms_28734_a
from primers.synthesis.test_models import primer_28727_left, primer_28727_right, \
    primer_28734_left, primer_28734_right
from targeted_enrichment.reagents.test_models import ter_28727, ter_28734

from ..runs.test_models import ngsrun
from lib_prep.workflows.test_models import magicalpcr1library, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a
from targeted_enrichment.unwrapping.test_models import pu_28727, pu_28734


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
