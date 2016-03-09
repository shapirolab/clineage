import pytest
import os
from Bio import SeqIO
import filecmp

from sequencing.analysis.adamiya import merge, create_reads_index, \
    align_primers_to_reads, _create_panel_fasta, _collect_mappings_from_sam, \
    _validate_unwrapper_mapping, _aggregate_read_ids_by_unwrapper, \
    seperate_reads_by_amplicons, _build_ms_variations, get_adam_ms_variations
from sequencing.analysis.models import AdamMSVariations, BowtieIndexMixin, \
    MicrosatelliteHistogramGenotype, name_to_ms_genotypes, ms_genotypes_to_name

index_files = ["{}.{}.bt2".format(BowtieIndexMixin.INDEX_PREFIX,x) for x in
    ["1","2","3","4","1.rev","2.rev"]]

def strip_fasta_records(fasta_records):
    for rec in fasta_records:
        yield rec.id, str(rec.seq)


@pytest.mark.django_db
def test_samplereads(samplereads):
    assert samplereads.barcoded_content.content.name == 'human amplified content'


@pytest.mark.django_db
def test_mergedreads(mergedreads, merged_reads_stripped_fasta):
    assert os.path.isfile(mergedreads.demux_reads.fastq1)
    assert os.path.isfile(mergedreads.demux_reads.fastq2)
    assert set(strip_fasta_records(mergedreads.included_reads_generator('M'))) == merged_reads_stripped_fasta
    assert set(strip_fasta_records(mergedreads.included_reads_generator('M'))) == set(strip_fasta_records(mergedreads.included_reads_generator('F')))


@pytest.mark.django_db
def test_runmerge(samplereads):
    mr = merge(samplereads)
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


@pytest.mark.django_db
def test_readsindex_bowtie2build(mergedreads):
    assert os.path.isfile(mergedreads.assembled_fastq)
    ri = create_reads_index(mergedreads, included_reads='M', padding=5)
    assert set(os.path.join(ri.index_dump_dir, x) for x in os.listdir(ri.index_dump_dir)) == set(ri.files)
    ri.delete()
    assert not os.path.exists(ri.index_dump_dir)


def test_create_panel_fasta(pu_28727, pu_28734):
    panel_fasta_name = _create_panel_fasta([pu_28727, pu_28734])
    assert os.path.isfile(panel_fasta_name)
    assert set(strip_fasta_records(
        SeqIO.parse(panel_fasta_name, "fasta")
    )) == {
        ("unwrapper_1_left", "TTTACTATGCCATGCTGCTGCT",),
        ("unwrapper_1_right", "TGTGCAAACAAGAACAGATGCC",),
        ("unwrapper_2_left", "AAGGCTTCTCCCCATTCCAAAG",),
        ("unwrapper_2_right", "AGTCCAAGCACACACTACTTCC",),
    }


def test_readsindex_unwrappers(readsindex_merged_only, pu_28727, pu_28734):
    assert list(
        readsindex_merged_only.merged_reads.demux_reads.library.unwrappers
    ) == [pu_28727, pu_28734]


@pytest.mark.django_db
def test_align_primers_to_reads_basic(readsindex_merged_only,
                                      require_unwrappers):
    ma = align_primers_to_reads(readsindex_merged_only)
    assert os.path.isfile(ma.assignment_sam)


@pytest.mark.django_db
def test_collect_mappings_from_sam(adammarginassignment,
                                   require_unwrappers,
                                   reads_matches):
    assert _collect_mappings_from_sam(adammarginassignment) == reads_matches


@pytest.mark.django_db
def test_align_primers_to_read_ids_with_mapping(readsindex_merged_only,
                                             adammarginassignment):
    assert \
        _collect_mappings_from_sam(
            align_primers_to_reads(readsindex_merged_only)
        ) == _collect_mappings_from_sam(adammarginassignment)


@pytest.mark.django_db
def test_validate_unwrapper_mapping(reads_matches, reads_unwrappers):
    assert set(_validate_unwrapper_mapping(reads_matches)) == reads_unwrappers


@pytest.mark.django_db
def test_aggregate_read_ids_by_unwrapper(reads_unwrappers, reads_by_unwrappers):
    assert _aggregate_read_ids_by_unwrapper(reads_unwrappers) == reads_by_unwrappers


@pytest.mark.django_db
def test_seperate_reads_by_amplicons(adammarginassignment, adamampliconreads):
    aars = list(seperate_reads_by_amplicons(adammarginassignment))
    assert len(aars) == 1
    aar = aars[0]
    aar_reads = set(strip_fasta_records(
        SeqIO.parse(aar.fastq, "fastq"))
    )
    ref_reads = set(strip_fasta_records(
        SeqIO.parse(adamampliconreads.fastq, "fastq"))
    )
    assert aar_reads == ref_reads

def test_build_ms_variations(pu_28727, pu_28727_adam_ms_variations):
    fasta = _build_ms_variations(pu_28727, 50)
    variations = set(strip_fasta_records(SeqIO.parse(fasta, "fasta")))
    assert variations == pu_28727_adam_ms_variations

@pytest.mark.django_db
def test_get_adam_ms_variations(pu_28727):
    assert AdamMSVariations.objects.count() == 0
    amsv = get_adam_ms_variations(pu_28727, 50)
    assert AdamMSVariations.objects.count() == 1
    amsv2 = get_adam_ms_variations(pu_28727, 50)
    assert AdamMSVariations.objects.count() == 1
    assert amsv2.id == amsv.id
    # TODO: test that we get a good index.
    amsv.delete()
    assert not os.path.exists(amsv.index_dump_dir)

@pytest.mark.django_db
def test_ms_histogram_genotype(ms_28727_a, ms_28727_b):
    assert MicrosatelliteHistogramGenotype.objects.count() == 0
    mhg1 = MicrosatelliteHistogramGenotype.get_for_genotype(ms_28727_a, 10)
    assert "{}".format(mhg1) == "1=10"
    assert MicrosatelliteHistogramGenotype.objects.count() == 1
    mhg2 = MicrosatelliteHistogramGenotype.get_for_string("1=10")
    assert MicrosatelliteHistogramGenotype.objects.count() == 1
    assert mhg1.id == mhg2.id
    assert mhg1.microsatellite == ms_28727_a
    assert mhg1.repeat_number == 10
    mhg3 = MicrosatelliteHistogramGenotype.get_for_string("2=20")
    assert MicrosatelliteHistogramGenotype.objects.count() == 2
    assert mhg3.microsatellite == ms_28727_b
    assert mhg3.repeat_number == 20
    assert mhg1.sequence == "TCT"*10
    assert mhg3.sequence == "CTG"*20
    mhg1.delete()
    mhg3.delete()


@pytest.mark.django_db
def test_ms_histogram_genotypes_names(ms_28727_a, ms_28734_a):
    mhg1 = MicrosatelliteHistogramGenotype.get_for_genotype(ms_28727_a, 5)
    mhg2 = MicrosatelliteHistogramGenotype.get_for_genotype(ms_28734_a, 8)
    assert ms_genotypes_to_name([mhg1,mhg2],prefix="abc") == "abc:1=5:3=8"
    mhgs, prefix = name_to_ms_genotypes("xyz:3=7:1=11")
    assert prefix == "xyz"
    mhg3, mhg4 = mhgs
    assert mhg3.microsatellite == ms_28734_a
    assert mhg3.repeat_number == 7
    assert mhg4.microsatellite == ms_28727_a
    assert mhg4.repeat_number == 11
    for mhg in [mhg1,mhg2,mhg3,mhg4]:
        mhg.delete()
