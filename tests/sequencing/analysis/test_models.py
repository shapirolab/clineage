import pytest
import os
from Bio import SeqIO

from sequencing.analysis.adamiya import merge, create_reads_index, \
    align_primers_to_reads, _create_panel_fasta, seperate_reads_by_amplicons, \
    get_adam_ms_variations, get_ms_variations_for_amplicon_reads, \
    separate_reads_by_genotypes
from sequencing.analysis.models import AdamMSVariations, \
    MicrosatelliteHistogramGenotype, name_to_ms_genotypes, ms_genotypes_to_name

from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM, \
    srs_to_tups, rc_srs_to_tups, strip_fasta_records
from tests.sequencing.analysis.reads_dict import ASSEMBLED, UNASSEMBLED


@pytest.mark.django_db
def test_samplereads(sample_reads_d):
    for sr in sample_reads_d.values():
        assert sr.barcoded_content.subclass.content.name == \
            'human amplified content'


@pytest.mark.django_db
def test_runmerge(sample_reads_d, adam_reads_fd):
    for (l_id, bc_id), sr in sample_reads_d.items():
        mr = merge(sr)
        assert mr.sample_reads is sr
        assert os.path.isfile(mr.assembled_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.assembled_fastq, "fastq"))) == \
            set(srs_to_tups(adam_reads_fd[l_id, bc_id, ASSEMBLED][RM]))
        assert os.path.isfile(mr.discarded_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.discarded_fastq, "fastq"))) == set()
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_forward_fastq, "fastq"))) == \
            set(srs_to_tups(adam_reads_fd[l_id, bc_id, UNASSEMBLED][R1]))
        assert os.path.isfile(mr.unassembled_reverse_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_reverse_fastq, "fastq"))) == \
            set(rc_srs_to_tups(adam_reads_fd[l_id, bc_id, UNASSEMBLED][R2]))
        mr.delete()
        assert not os.path.exists(mr.assembled_fastq)
        assert not os.path.exists(mr.discarded_fastq)
        assert not os.path.exists(mr.unassembled_forward_fastq)
        assert not os.path.exists(mr.unassembled_reverse_fastq)


@pytest.mark.django_db
def test_merged_reads(adam_merged_reads_d, adam_reads_fd):
    for (l_id, bc_id), mr in adam_merged_reads_d.items():
        assert set(srs_to_tups(mr.included_reads_generator('M'))) == \
             set(srs_to_tups(adam_reads_fd[l_id, bc_id, ASSEMBLED][RM]))
        assert set(srs_to_tups(mr.included_reads_generator('F'))) == \
             set(srs_to_tups(adam_reads_fd[l_id, bc_id, ASSEMBLED][RM])) | \
             set(srs_to_tups(adam_reads_fd[l_id, bc_id, UNASSEMBLED][R1]))


def test_create_panel_fasta(pu_28727, pu_28734):
    panel_fasta_name = _create_panel_fasta([pu_28727, pu_28734])
    assert os.path.isfile(panel_fasta_name)
    assert set(strip_fasta_records(
        SeqIO.parse(panel_fasta_name, "fasta")
    )) == {
        ("amplicon_28727_left", "TTTACTATGCCATGCTGCTGCT",),
        ("amplicon_28727_right", "TGTGCAAACAAGAACAGATGCC",),
        ("amplicon_28734_left", "AAGGCTTCTCCCCATTCCAAAG",),
        ("amplicon_28734_right", "AGTCCAAGCACACACTACTTCC",),
    }
    os.unlink(panel_fasta_name)


def test_amplicons_mapping(adam_merged_reads_d, adam_reads_fd, requires_amplicons):
    for (l_id, bc_id), mr in adam_merged_reads_d.items():
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
                amp = aar.amplicon_id
                amps.add(amp)
                aar_fnames_d = {
                    R1: aar.fastq1,
                    R2: aar.fastq2,
                    RM: aar.fastqm,
                }
                if inc == "M":
                    ref_reads_d = adam_reads_fd[l_id, bc_id, ASSEMBLED, amp]
                else:  # inc == "F"
                    ref_reads_d = {
                        R1: adam_reads_fd[l_id, bc_id, ASSEMBLED, amp][R1] + \
                            adam_reads_fd[l_id, bc_id, UNASSEMBLED, amp][R1],
                        R2: adam_reads_fd[l_id, bc_id, ASSEMBLED, amp][R2] + \
                            adam_reads_fd[l_id, bc_id, UNASSEMBLED, amp][R2],
                        RM: adam_reads_fd[l_id, bc_id, ASSEMBLED, amp][RM] + \
                            adam_reads_fd[l_id, bc_id, UNASSEMBLED, amp][R1],
                    }
                for r in [R1, R2, RM]:
                    assert set(srs_to_tups(
                        SeqIO.parse(aar_fnames_d[r], "fastq"))
                    ) == \
                    set(srs_to_tups(
                        ref_reads_d[r]
                    ))
                aar.delete()
                assert not os.path.exists(aar.fastq1)
                assert not os.path.exists(aar.fastq2)
                assert not os.path.exists(aar.fastqm)
            assert amps == set(adam_reads_fd.sub((l_id, bc_id, ASSEMBLED)).keys()) | \
                (set(adam_reads_fd.sub((l_id, bc_id, UNASSEMBLED)).keys()) if \
                    inc == "F" else set())
            ama.delete()
            assert not os.path.exists(ama.assignment_sam)
            ri.delete()
            assert not os.path.exists(ri.index_dump_dir)


def test_genotype_mapping(adam_amplicon_reads_d, adam_reads_fd, requires_amplicons, requires_microsatellites):
    for (l_id, bc, inc, amp), aar in adam_amplicon_reads_d.items():
        ms_planning_version = 0
        padding = 50
        # FIXME
        # test_align_reads_to_ms_variations
        ah = get_ms_variations_for_amplicon_reads(aar, padding, ms_planning_version)
        assert ah.amplicon_reads_id == aar.id
        # TODO: make this check more explicit.
        assert ah.ms_variations == AdamMSVariations.objects.get(
            amplicon_id=amp,
            padding=padding,
            microsatellites_version=ms_planning_version
        )
        assert ah.amplicon_id == amp
        assert ah.microsatellites_version == ms_planning_version
        assert os.path.isfile(ah.assignment_sam)
        gens = set()
        # test_separate_reads_by_genotypes
        for her in separate_reads_by_genotypes(ah):
            assert her.histogram_id == ah.id
            assert set(her.snp_genotypes.genotypes) == set()
            gen = frozenset((msg.microsatellite_id, msg.repeat_number) for \
                msg in her.microsatellite_genotypes.genotypes)
            gens.add(gen)

            her_fnames_d = {
                R1: her.fastq1,
                R2: her.fastq2,
                RM: her.fastqm,
            }
            if inc == "M":
                ref_reads_d = adam_reads_fd[l_id, bc, ASSEMBLED, amp, gen]
            else:  # inc == "F"
                ref_reads_d = {
                    R1: adam_reads_fd[l_id, bc, ASSEMBLED, amp, gen][R1] + \
                        adam_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R1],
                    R2: adam_reads_fd[l_id, bc, ASSEMBLED, amp, gen][R2] + \
                        adam_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R2],
                    RM: adam_reads_fd[l_id, bc, ASSEMBLED, amp, gen][RM] + \
                        adam_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R1],
                }
            for r in [R1, R2, RM]:
                assert set(srs_to_tups(  # TODO: get informative error on genotyping mismatch
                    SeqIO.parse(her_fnames_d[r], "fastq"))
                ) == \
                set(srs_to_tups(
                    ref_reads_d[r]
                ))
            assert her.num_reads == \
                len(ref_reads_d[RM])
            her.delete()
            assert not os.path.exists(her.fastq1)
            assert not os.path.exists(her.fastq2)
            assert not os.path.exists(her.fastqm)
        if inc == "F":
            assert gens == set(adam_reads_fd.sub((l_id, bc, ASSEMBLED, amp)).keys()) \
                | set(adam_reads_fd.sub((l_id, bc, UNASSEMBLED, amp)).keys())
        else:
            assert gens == set(adam_reads_fd.sub((l_id, bc, ASSEMBLED, amp)).keys())
        ah.delete()
        assert not os.path.exists(ah.assignment_sam)
        #assert filecmp.cmp(ah.assignment_sam,
    # FIXME: get these from outside.
    AdamMSVariations.objects.all().delete()


def test_get_adam_ms_variations(pu_28727, pu_28734, requires_microsatellites):
    ms_planning_version = 0
    assert AdamMSVariations.objects.count() == 0
    amsv1 = get_adam_ms_variations(pu_28727, 50, ms_planning_version)
    assert amsv1.padding == 50
    assert amsv1.amplicon_id == pu_28727.id
    assert amsv1.microsatellites_version == ms_planning_version
    assert os.path.isdir(amsv1.index_dump_dir)
    assert AdamMSVariations.objects.count() == 1
    amsv2 = get_adam_ms_variations(pu_28727, 30, ms_planning_version)
    assert amsv2.padding == 30
    assert amsv2.amplicon_id == pu_28727.id
    assert amsv2.microsatellites_version == ms_planning_version
    assert os.path.isdir(amsv2.index_dump_dir)
    assert AdamMSVariations.objects.count() == 2
    amsv3 = get_adam_ms_variations(pu_28727, 50, ms_planning_version)
    assert amsv3.padding == 50
    assert amsv3.amplicon_id == pu_28727.id
    assert amsv3.microsatellites_version == ms_planning_version
    assert os.path.isdir(amsv3.index_dump_dir)
    assert AdamMSVariations.objects.count() == 2
    amsv4 = get_adam_ms_variations(pu_28734, 50, ms_planning_version)
    assert amsv4.padding == 50
    assert amsv4.amplicon_id == pu_28734.id
    assert amsv4.microsatellites_version == ms_planning_version
    assert os.path.isdir(amsv4.index_dump_dir)
    assert AdamMSVariations.objects.count() == 3
    assert amsv1.id != amsv2.id
    assert amsv1.id != amsv4.id
    assert amsv2.id != amsv4.id
    assert amsv1.id == amsv3.id
    amsv3.delete()
    assert not os.path.exists(amsv3.index_dump_dir)
    assert AdamMSVariations.objects.count() == 2
    amsv5 = get_adam_ms_variations(pu_28727, 50, ms_planning_version)
    assert AdamMSVariations.objects.count() == 3
    assert amsv5.padding == 50
    assert amsv5.amplicon_id == pu_28727.id
    assert amsv5.microsatellites_version == ms_planning_version
    assert os.path.isdir(amsv5.index_dump_dir)
    amsv2.delete()
    assert not os.path.exists(amsv2.index_dump_dir)
    amsv4.delete()
    assert not os.path.exists(amsv4.index_dump_dir)
    amsv5.delete()
    assert not os.path.exists(amsv5.index_dump_dir)


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
