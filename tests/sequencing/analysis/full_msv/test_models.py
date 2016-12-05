import pytest
import itertools
import os
import io
import plumbum
from Bio import SeqIO

from sequencing.analysis.full_msv.full_msv import merge, \
    get_full_ms_variations, separate_reads_by_genotypes, \
    align_reads_to_ms_variations, stream_group_alignemnts
from sequencing.analysis.models import FullMSVHistogram, HistogramEntryReads

from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM, \
    srs_to_tups, rc_srs_to_tups, strip_fasta_records
from tests.sequencing.analysis.full_msv.reads_dict import ASSEMBLED, UNASSEMBLED


@pytest.mark.django_db
def test_samplereads(sample_reads_d):
    for sr in sample_reads_d.values():
        assert sr.barcoded_content.subclass.content.name == \
            'human amplified content'


@pytest.mark.django_db
def test_runmerge(sample_reads_d, fmsv_reads_fd):
    for (l_id, bc_id), sr in sample_reads_d.items():
        mr = merge(sr)
        assert mr.sample_reads is sr
        assert os.path.isfile(mr.assembled_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.assembled_fastq, "fastq"))) == \
            set(srs_to_tups(fmsv_reads_fd[l_id, bc_id, ASSEMBLED][RM]))
        assert os.path.isfile(mr.discarded_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.discarded_fastq, "fastq"))) == set()
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_forward_fastq, "fastq"))) == \
            set(srs_to_tups(fmsv_reads_fd[l_id, bc_id, UNASSEMBLED][R1]))
        assert os.path.isfile(mr.unassembled_reverse_fastq)
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_reverse_fastq, "fastq"))) == \
            set(rc_srs_to_tups(fmsv_reads_fd[l_id, bc_id, UNASSEMBLED][R2]))
        mr.delete()
        assert not os.path.exists(mr.assembled_fastq)
        assert not os.path.exists(mr.discarded_fastq)
        assert not os.path.exists(mr.unassembled_forward_fastq)
        assert not os.path.exists(mr.unassembled_reverse_fastq)


@pytest.mark.django_db
def test_merged_reads(fmsv_merged_reads_d, fmsv_reads_fd):
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        assert set(srs_to_tups(mr.included_reads_generator('M'))) == \
             set(srs_to_tups(fmsv_reads_fd[l_id, bc_id, ASSEMBLED][RM]))
        assert set(srs_to_tups(mr.included_reads_generator('F'))) == \
            set(srs_to_tups(fmsv_reads_fd[l_id, bc_id, ASSEMBLED][RM])) | \
            set(srs_to_tups(fmsv_reads_fd[l_id, bc_id, UNASSEMBLED][R1]))


@pytest.mark.django_db
def test_build_ms_variations(amp_134_full_msv, ac_134, requires_microsatellites):
    padding = 50
    mss_version = 0
    amplicon_collection = ac_134
    fmsv = get_full_ms_variations(amplicon_collection, padding, mss_version)
    bowtie2_inspect = plumbum.local["bowtie2-inspect"]
    s = bowtie2_inspect(fmsv.index_files_prefix)
    fasta = io.StringIO(s)
    variations = set(strip_fasta_records(SeqIO.parse(fasta, "fasta")))
    assert variations == amp_134_full_msv
    fasta.close()
    fmsv.delete()
    assert not os.path.exists(fmsv.index_dump_dir)


@pytest.mark.django_db
def test_stream_group_alignemnts(fmsv_merged_reads_d, fmsv_reads_fd, full_ms_variations, amplicon_collection, requires_none_genotypes):
    padding = 50
    mss_version = 0
    msv = get_full_ms_variations(amplicon_collection, padding, mss_version)
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        # for inc in ["M", "F"]:
        for inc in ["M"]:
            # test_ms_variations_index
            assert os.path.isfile(mr.assembled_fastq)

            fmsva = align_reads_to_ms_variations(mr, padding, mss_version)  #included_reads=inc
            for amp_id, histogram_reads in stream_group_alignemnts(fmsva):
                generator_read_ids = []
                for read_ids in histogram_reads.values():
                    generator_read_ids += read_ids
                if inc == "M":
                    ref_reads_d = [seq.id for seq in fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp_id][RM]]
                # TODO: inc == "F"
                assert set(generator_read_ids) == set(ref_reads_d)
            fmsva.delete()
    msv.delete()


@pytest.mark.django_db
def test_amplicons_mapping(fmsv_merged_reads_d, fmsv_reads_fd, full_ms_variations, amplicon_collection, requires_none_genotypes):
    padding = 50
    mss_version = 0
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        # for inc in ["M", "F"]:
        for inc in ["M"]:
            # test_ms_variations_index
            assert os.path.isfile(mr.assembled_fastq)
            fmsva = align_reads_to_ms_variations(mr, padding, mss_version)  #included_reads=inc
            assert fmsva.merged_reads == mr
            assert set(os.path.join(fmsva.ms_variations.index_dump_dir, x) for x in \
                       os.listdir(fmsva.ms_variations.index_dump_dir)) == set(fmsva.ms_variations.files)

            # test_separate_reads_by_amplicons
            for her in separate_reads_by_genotypes(fmsva):
                pass
            amps = set()
            assert FullMSVHistogram.objects.filter(
                assignment=fmsva,
            ).count() > 0
            for fmsvh in FullMSVHistogram.objects.filter(
                assignment=fmsva,
            ):
                fmsvh_reads = {
                    'R1': itertools.chain(*[
                        SeqIO.parse(her.fastq1, "fastq") for her
                        in HistogramEntryReads.objects.filter(histogram=fmsvh)
                    ]),
                    'R2': itertools.chain(*[
                        SeqIO.parse(her.fastq2, "fastq") for her
                        in HistogramEntryReads.objects.filter(histogram=fmsvh)
                    ]),
                    'RM': itertools.chain(*[
                        SeqIO.parse(her.fastqm, "fastq") for her
                        in HistogramEntryReads.objects.filter(histogram=fmsvh)
                    ])
                }
                amp = fmsvh.amplicon_id
                amps.add(amp)
                if inc == "M":
                    ref_reads_d = fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp]
                else:  # inc == "F"
                    ref_reads_d = {
                        R1: fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp][R1] +
                            fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp][R1],
                        R2: fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp][R2] +
                            fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp][R2],
                        RM: fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp][RM] +
                            fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp][R1],
                    }
                for r in [R1, R2, RM]:
                    assert set(srs_to_tups(
                        fmsvh_reads[r]
                    )) == set(srs_to_tups(
                        ref_reads_d[r]
                    ))
                for her in HistogramEntryReads.objects.filter(histogram=fmsvh):
                    her.delete()
                    assert not os.path.exists(her.fastq1)
                    assert not os.path.exists(her.fastq2)
                    assert not os.path.exists(her.fastqm)
            assert amps == set(fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED)) | \
                           (set(fmsv_reads_fd.keys(l_id, bc_id, UNASSEMBLED)) if \
                                inc == "F" else set())
            fmsvv = fmsva.ms_variations
            fmsva.delete()
            assert not os.path.exists(fmsva.sorted_assignment_bam)
            fmsvv.delete()
            assert not os.path.exists(fmsvv.index_dump_dir)


@pytest.mark.django_db
def test_genotype_mapping(fmsv_merged_reads_d, fmsv_reads_fd, full_ms_variations, amplicon_collection, requires_none_genotypes):
    padding = 50
    mss_version = 0
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        # for inc in ["M", "F"]:
        for inc in ["M"]:
            # test_ms_variations_index
            assert os.path.isfile(mr.assembled_fastq)
            fmsva = align_reads_to_ms_variations(mr, padding, mss_version)  #included_reads=inc
            assert fmsva.merged_reads == mr
            assert set(os.path.join(fmsva.ms_variations.index_dump_dir, x) for x in \
                       os.listdir(fmsva.ms_variations.index_dump_dir)) == set(fmsva.ms_variations.files)
            amps = set()
            gens = set()
            # test_separate_reads_by_amplicons
            for her in separate_reads_by_genotypes(fmsva):
                # assert her.histogram_id == ah.id
                amp = her.histogram.amplicon_id
                amps.add(amp)
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
                    ref_reads_d = fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen]
                else:  # inc == "F"
                    ref_reads_d = {
                        R1: fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen][R1] + \
                            fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp, gen][R1],
                        R2: fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen][R2] + \
                            fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp, gen][R2],
                        RM: fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen][RM] + \
                            fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp, gen][R1],
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
            assert gens == set(fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED, amp)) \
                           | set(fmsv_reads_fd.keys(l_id, bc_id, UNASSEMBLED, amp))
        else:
            ref_gens = set()
            for amp in fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED):
                for s in fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED, amp):
                    ref_gens.add(s)
            assert gens == ref_gens
            assert amps == set(fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED)) | \
                           (set(fmsv_reads_fd.keys(l_id, bc_id, UNASSEMBLED)) if \
                                inc == "F" else set())
        fmsvv = fmsva.ms_variations
        fmsva.delete()
        assert not os.path.exists(fmsva.sorted_assignment_bam)
        fmsvv.delete()
        assert not os.path.exists(fmsvv.index_dump_dir)
