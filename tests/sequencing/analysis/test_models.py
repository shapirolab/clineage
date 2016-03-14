import pytest
import os
from Bio import SeqIO
import filecmp

from sequencing.analysis.adamiya import merge, create_reads_index, \
    align_primers_to_reads, _create_panel_fasta, _collect_mappings_from_sam, \
    _validate_amplicon_mapping, _aggregate_read_ids_by_amplicon, \
    seperate_reads_by_amplicons, _build_ms_variations, \
    get_adam_ms_variations, align_reads_to_ms_variations, \
    separate_reads_by_genotypes
from sequencing.analysis.models import AdamMSVariations, BowtieIndexMixin, \
    MicrosatelliteHistogramGenotype, name_to_ms_genotypes, ms_genotypes_to_name

index_files = ["{}.{}.bt2".format(BowtieIndexMixin.INDEX_PREFIX,x) for x in
    ["1","2","3","4","1.rev","2.rev"]]

def strip_fasta_records(fasta_records):
    for rec in fasta_records:
        yield rec.id, str(rec.seq)


@pytest.mark.django_db
def test_samplereads(samplereads_bc1):
    assert samplereads_bc1.barcoded_content.content.name == 'human amplified content'


@pytest.mark.django_db
def test_mergedreads(mergedreads_bc1, merged_reads_stripped_fasta):
    assert os.path.isfile(mergedreads_bc1.sample_reads.fastq1)
    assert os.path.isfile(mergedreads_bc1.sample_reads.fastq2)
    assert set(strip_fasta_records(mergedreads_bc1.included_reads_generator('M'))) == merged_reads_stripped_fasta
    assert set(strip_fasta_records(mergedreads_bc1.included_reads_generator('M'))) == set(strip_fasta_records(mergedreads_bc1.included_reads_generator('F')))


@pytest.mark.django_db
def test_runmerge(samplereads_bc1):
    mr = merge(samplereads_bc1)
    assert mr.sample_reads is samplereads_bc1
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
def test_readsindex_bowtie2build(mergedreads_bc1):
    assert os.path.isfile(mergedreads_bc1.assembled_fastq)
    ri = create_reads_index(mergedreads_bc1, included_reads='M', padding=5)
    assert set(os.path.join(ri.index_dump_dir, x) for x in os.listdir(ri.index_dump_dir)) == set(ri.files)
    ri.delete()
    assert not os.path.exists(ri.index_dump_dir)


def test_create_panel_fasta(pu_28727, pu_28734):
    panel_fasta_name = _create_panel_fasta([pu_28727, pu_28734])
    assert os.path.isfile(panel_fasta_name)
    assert set(strip_fasta_records(
        SeqIO.parse(panel_fasta_name, "fasta")
    )) == {
        ("amplicon_1_left", "TTTACTATGCCATGCTGCTGCT",),
        ("amplicon_1_right", "TGTGCAAACAAGAACAGATGCC",),
        ("amplicon_2_left", "AAGGCTTCTCCCCATTCCAAAG",),
        ("amplicon_2_right", "AGTCCAAGCACACACTACTTCC",),
    }


def test_readsindex_amplicons(readsindex_bc1_M, pu_28727, pu_28734):
    assert list(
        readsindex_bc1_M.merged_reads.sample_reads.library.amplicons
    ) == [pu_28727, pu_28734]


@pytest.mark.django_db
def test_align_primers_to_reads_basic(readsindex_bc1_M,
                                      require_amplicons):
    ma = align_primers_to_reads(readsindex_bc1_M)
    assert os.path.isfile(ma.assignment_sam)


@pytest.mark.django_db
def test_collect_mappings_from_sam(adammarginassignment_bc1_F,
                                   require_amplicons,
                                   reads_matches):
    assert _collect_mappings_from_sam(adammarginassignment_bc1_F) == reads_matches


@pytest.mark.django_db
def test_align_primers_to_read_ids_with_mapping(readsindex_bc1_M,
                                             adammarginassignment_bc1_F):
    assert \
        _collect_mappings_from_sam(
            align_primers_to_reads(readsindex_bc1_M)
        ) == _collect_mappings_from_sam(adammarginassignment_bc1_F)


@pytest.mark.django_db
def test_validate_amplicon_mapping(reads_matches, reads_amplicons):
    assert set(_validate_amplicon_mapping(reads_matches)) == reads_amplicons


@pytest.mark.django_db
def test_aggregate_read_ids_by_amplicon(reads_amplicons, reads_by_amplicons):
    assert _aggregate_read_ids_by_amplicon(reads_amplicons) == reads_by_amplicons


@pytest.mark.django_db
def test_seperate_reads_by_amplicons(adammarginassignment_bc1_F, adamampliconreads_bc1_F_28734):
    aars = list(seperate_reads_by_amplicons(adammarginassignment_bc1_F))
    assert len(aars) == 1
    aar = aars[0]
    aar_reads = set(strip_fasta_records(
        SeqIO.parse(aar.fastq, "fastq"))
    )
    ref_reads = set(strip_fasta_records(
        SeqIO.parse(adamampliconreads_bc1_F_28734.fastq, "fastq"))
    )
    assert aar_reads == ref_reads
    aar_reads1 = set(strip_fasta_records(
        SeqIO.parse(aar.fastq1, "fastq"))
    )
    ref_reads1 = set(strip_fasta_records(
        SeqIO.parse(adamampliconreads_bc1_F_28734.fastq1, "fastq"))
    )
    assert aar_reads1 == ref_reads1
    aar_reads2 = set(strip_fasta_records(
        SeqIO.parse(aar.fastq2, "fastq"))
    )
    ref_reads2 = set(strip_fasta_records(
        SeqIO.parse(adamampliconreads_bc1_F_28734.fastq2, "fastq"))
    )
    assert aar_reads2 == ref_reads2

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

@pytest.mark.django_db
def test_align_reads_to_ms_variations(adamampliconreads_bc1_F_28734):
    ah = align_reads_to_ms_variations(adamampliconreads_bc1_F_28734, 50)
    assert os.path.isfile(ah.assignment_sam)
    #assert filecmp.cmp(ah.assignment_sam, 

@pytest.mark.django_db
def test_separate_reads_by_genotypes(adamhistogram, pu_28734, ms_28734_a):
    l = list(separate_reads_by_genotypes(adamhistogram))
    assert len(l) == 1
    her = l[0]
    assert her.histogram == adamhistogram
    assert her.amplicon == pu_28734
    assert set(her.microsatellite_genotypes.all()) == \
        {MicrosatelliteHistogramGenotype.get_for_genotype(ms_28734_a, 7)}
    assert set(her.snp_genotypes.all()) == set()
    assert her.num_reads == 12
    #fastq1 = models.FilePathField()
    #fastq2 = models.FilePathField()
    #fastqm = models.FilePathField(null=True)
