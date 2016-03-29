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

from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM, \
    srs_to_tups, rc_srs_to_tups
from tests.sequencing.analysis.reads_dict import ASSEMBLED, UNASSEMBLED

index_files = ["{}.{}.bt2".format(BowtieIndexMixin.INDEX_PREFIX,x) for x in
    ["1","2","3","4","1.rev","2.rev"]]

def strip_fasta_records(fasta_records):
    for rec in fasta_records:
        yield rec.id, str(rec.seq)


@pytest.mark.django_db
def test_samplereads(sample_reads_d):
    for sr in sample_reads_d.itervalues():
        assert sr.barcoded_content.content.name == 'human amplified content'


@pytest.mark.django_db
def test_runmerge(sample_reads_d, adam_reads_fd):
    for k, sr in sample_reads_d.iteritems():
        mr = merge(sr)
        assert mr.sample_reads is sr
        assert os.path.isfile(mr.assembled_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.assembled_fastq, "fastq"))) == \
            set(srs_to_tups(adam_reads_fd[k,ASSEMBLED][RM]))
        assert os.path.isfile(mr.discarded_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.discarded_fastq, "fastq"))) == set()
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_forward_fastq, "fastq"))) == \
            set(srs_to_tups(adam_reads_fd[k,UNASSEMBLED][R1]))
        assert os.path.isfile(mr.unassembled_reverse_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_reverse_fastq, "fastq"))) == \
            set(rc_srs_to_tups(adam_reads_fd[k,UNASSEMBLED][R2]))
        mr.delete()
        assert not os.path.exists(mr.assembled_fastq)
        assert not os.path.exists(mr.discarded_fastq)
        assert not os.path.exists(mr.unassembled_forward_fastq)
        assert not os.path.exists(mr.unassembled_reverse_fastq)


@pytest.mark.django_db
def test_merged_reads(adam_merged_reads_d, adam_reads_fd):
    for k, mr in adam_merged_reads_d.iteritems():
        assert set(srs_to_tups(mr.included_reads_generator('M'))) == \
             set(srs_to_tups(adam_reads_fd[k,ASSEMBLED][RM]))
        assert set(srs_to_tups(mr.included_reads_generator('F'))) == \
             set(srs_to_tups(adam_reads_fd[k,ASSEMBLED][RM])) | \
             set(srs_to_tups(adam_reads_fd[k,UNASSEMBLED][R1]))


#@pytest.mark.django_db
#def test_runmerge_custom_read(samplereads_bc2):
    #mr = merge(samplereads_bc2)
    #reads = list(strip_fasta_records(mr.included_reads_generator('M')))
    #assert len(reads) == 1
    #read = reads[0]
    #assert read == (
            #'M00393:167:000000000-ABF3N:1:1101:20583:16769',
            #'AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        #)
    #mr.delete()

#@pytest.mark.django_db
#def test_adamhistogram_plus3(adamhistogram_plus3, pu_28734, ms_28734_a):
    #l = list(separate_reads_by_genotypes(adamhistogram_plus3))
    #assert len(l) == 1
    #her = l[0]
    #assert her.histogram == adamhistogram_plus3
    #assert her.amplicon == pu_28734
    #assert set(her.microsatellite_genotypes.all()) == \
        #{MicrosatelliteHistogramGenotype.get_for_genotype(ms_28734_a, 9)}
    #assert set(her.snp_genotypes.all()) == set()
    #assert her.num_reads == 1


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
    os.unlink(panel_fasta_name)

def test_amplicons_mapping(adam_merged_reads_d, adam_reads_fd, amplicon_d_r):
    for bc, mr in adam_merged_reads_d.iteritems():
        for inc in ["M", "F"]:
            # test_readsindex_bowtie2build
            assert os.path.isfile(mr.assembled_fastq)
            ri = create_reads_index(mr, included_reads=inc, padding=5)
            assert ri.merged_reads == mr
            assert set(os.path.join(ri.index_dump_dir, x) for x in \
                os.listdir(ri.index_dump_dir)) == set(ri.files)

            #test_align_primers_to_reads_basic
            ama = align_primers_to_reads(ri)
            assert os.path.isfile(ama.assignment_sam)

            #test_seperate_reads_by_amplicons
            amps = set()
            for aar in seperate_reads_by_amplicons(ama):
                amp = amplicon_d_r[aar.amplicon]
                amps.add(amp)
                aar_fnames_d = {
                    R1: aar.fastq1,
                    R2: aar.fastq2,
                    RM: aar.fastqm,
                }
                if inc == "M":
                    ref_reads_d = adam_reads_fd[bc, ASSEMBLED, amp]
                else:  # inc == "F"
                    ref_reads_d = {
                        R1: adam_reads_fd[bc, ASSEMBLED, amp][R1] + \
                            adam_reads_fd[bc, UNASSEMBLED, amp][R1],
                        R2: adam_reads_fd[bc, ASSEMBLED, amp][R2] + \
                            adam_reads_fd[bc, UNASSEMBLED, amp][R2],
                        RM: adam_reads_fd[bc, ASSEMBLED, amp][RM] + \
                            adam_reads_fd[bc, UNASSEMBLED, amp][R1],
                    }
                for r in [R1, R2, RM]:
                    assert set(strip_fasta_records(
                        SeqIO.parse(aar_fnames_d[r], "fastq"))
                    ) == \
                    set(strip_fasta_records(
                        ref_reads_d[r]
                    ))
                aar.delete()
                assert not os.path.exists(aar.fastq1)
                assert not os.path.exists(aar.fastq2)
                assert not os.path.exists(aar.fastqm)
            assert amps == set(adam_reads_fd.sub((bc, ASSEMBLED)).keys())
            ama.delete()
            assert not os.path.exists(ama.assignment_sam)
            ri.delete()
            assert not os.path.exists(ri.index_dump_dir)


def test_genotype_mapping(adam_amplicon_reads_d, adam_reads_fd, amplicon_d_r):
    for (bc, inc, amp), amr in adam_amplicon_reads_d.iteritems():
        # FIXME
        # test_align_reads_to_ms_variations
        ah = align_reads_to_ms_variations(amr, 50)
        assert ah.amplicon_reads == amr
        assert os.path.isfile(ah.assignment_sam)
        gens = set()
        # test_separate_reads_by_genotypes
        for her in separate_reads_by_genotypes(ah):
            assert her.histogram_id == ah.id
            assert amplicon_d_r[her.amplicon] == amp
            assert set(her.snp_genotypes.all()) == set()
            gen = frozenset((msg.microsatellite_id, msg.repeat_number) for \
                msg in her.microsatellite_genotypes.all())
            gens.add(gen)

            her_fnames_d = {
                R1: her.fastq1,
                R2: her.fastq2,
                RM: her.fastqm,
            }
            if inc == "M":
                ref_reads_d = adam_reads_fd[bc, ASSEMBLED, amp, gen]
            else:  # inc == "F"
                ref_reads_d = {
                    R1: adam_reads_fd[bc, ASSEMBLED, amp, gen][R1] + \
                        adam_reads_fd[bc, UNASSEMBLED, amp, gen][R1],
                    R2: adam_reads_fd[bc, ASSEMBLED, amp, gen][R2] + \
                        adam_reads_fd[bc, UNASSEMBLED, amp, gen][R2],
                    RM: adam_reads_fd[bc, ASSEMBLED, amp, gen][RM] + \
                        adam_reads_fd[bc, UNASSEMBLED, amp, gen][R1],
                }
            for r in [R1, R2, RM]:
                assert set(strip_fasta_records(
                    SeqIO.parse(her_fnames_d[r], "fastq"))
                ) == \
                set(strip_fasta_records(
                    ref_reads_d[r]
                ))
            assert her.num_reads == \
                len(ref_reads_d[RM])
            her.delete()
            assert not os.path.exists(her.fastq1)
            assert not os.path.exists(her.fastq2)
            assert not os.path.exists(her.fastqm)
        ah.delete()
        assert not os.path.exists(ah.assignment_sam)
        #assert filecmp.cmp(ah.assignment_sam,
    # FIXME: get these from outside.
    AdamMSVariations.objects.all().delete()


def test_get_adam_ms_variations(pu_28727, pu_28734):
    assert AdamMSVariations.objects.count() == 0
    amsv1 = get_adam_ms_variations(pu_28727, 50)
    assert amsv1.padding == 50
    assert amsv1.amplicon == pu_28727
    assert os.path.isdir(amsv1.index_dump_dir)
    assert AdamMSVariations.objects.count() == 1
    amsv2 = get_adam_ms_variations(pu_28727, 30)
    assert amsv2.padding == 30
    assert amsv2.amplicon == pu_28727
    assert os.path.isdir(amsv2.index_dump_dir)
    assert AdamMSVariations.objects.count() == 2
    amsv3 = get_adam_ms_variations(pu_28727, 50)
    assert amsv3.padding == 50
    assert amsv3.amplicon == pu_28727
    assert os.path.isdir(amsv3.index_dump_dir)
    assert AdamMSVariations.objects.count() == 2
    amsv4 = get_adam_ms_variations(pu_28734, 50)
    assert amsv4.padding == 50
    assert amsv4.amplicon == pu_28734
    assert os.path.isdir(amsv4.index_dump_dir)
    assert AdamMSVariations.objects.count() == 3
    assert amsv1.id != amsv2.id
    assert amsv1.id != amsv4.id
    assert amsv2.id != amsv4.id
    assert amsv1.id == amsv3.id
    amsv3.delete()
    assert not os.path.exists(amsv3.index_dump_dir)
    assert AdamMSVariations.objects.count() == 2
    amsv5 = get_adam_ms_variations(pu_28727, 50)
    assert AdamMSVariations.objects.count() == 3
    assert amsv5.padding == 50
    assert amsv5.amplicon == pu_28727
    assert os.path.isdir(amsv5.index_dump_dir)
    amsv2.delete()
    assert not os.path.exists(amsv2.index_dump_dir)
    amsv4.delete()
    assert not os.path.exists(amsv4.index_dump_dir)
    amsv5.delete()
    assert not os.path.exists(amsv5.index_dump_dir)


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
