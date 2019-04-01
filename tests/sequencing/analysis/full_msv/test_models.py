import pytest
import itertools
import os
import io
import plumbum
from Bio import SeqIO

from utils.groups import grouper
from targeted_enrichment.amplicons.models import AmpliconCollection
from sequencing.analysis.full_msv.full_msv import merge, \
    get_full_ms_variations, separate_reads_by_genotypes, \
    align_reads_to_ms_variations, stream_group_alignemnts, \
    split_merged_reads, align_reads_to_ms_variations_part, \
    merge_fmsva_parts, split_merged_reads_as_list, align_reads_to_ms_variations_part,\
    merge_fmsva_parts_as_list, align_reads_to_ms_variations_part_as_list
from sequencing.analysis.models import HistogramEntryReads
from sequencing.analysis.full_msv.models import FullMSVHistogram,  \
        FullMSVariations, FullMSVAssignment, FullMSVAssignmentPart, FullMSVMergedReads, \
        FullMSVMergedReadsPart

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
def test_split_merged_reads(fmsv_merged_reads_d, fmsv_reads_fd):
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        for mrp in split_merged_reads(mr, 1, included_reads='M'):
            assert mrp.merged_reads == mr
        joined_reads = []
        for mrp in mr.fullmsvmergedreadspart_set.all():
            joined_reads += list(SeqIO.parse(mrp.fastq_part, "fastq"))
        assert set(srs_to_tups(mr.included_reads_generator('M'))) == set(srs_to_tups(joined_reads))
        mr.fullmsvmergedreadspart_set.all().delete()


@pytest.mark.django_db
def test_align_split_merged_reads(fmsv_merged_reads_d, fmsv_reads_fd, ac_134, requires_microsatellites):
    padding = 50
    mss_version = 0
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        amplicon_collection = mr.sample_reads.library.subclass.panel.amplicon_collection
        fmsv = get_full_ms_variations(amplicon_collection, padding, mss_version)
        for mrp in split_merged_reads(mr, 1, included_reads='M'):
            fmsvaps = align_reads_to_ms_variations_part_as_list(mrp, padding, mss_version)
            for fmsvap in fmsvaps:
                assert fmsvap.merged_reads_part == mrp
                assert os.path.isfile(fmsvap.assignment_bam)
        mr.fullmsvmergedreadspart_set.all().delete()
        fmsv.delete()


@pytest.mark.django_db
def test_colliding_split_merged_reads(fmsv_merged_reads_d, fmsv_reads_fd, ac_134, requires_microsatellites):
    padding = 50
    mss_version = 0
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        amplicon_collection = mr.sample_reads.library.subclass.panel.amplicon_collection
        fmsv = get_full_ms_variations(amplicon_collection, padding, mss_version)
        for mrp in split_merged_reads(mr, 1, included_reads='M'):
            fmsvaps = align_reads_to_ms_variations_part_as_list(mrp, padding, mss_version)
            for fmsvap in fmsvaps:
                assert fmsvap.merged_reads_part == mrp
                assert os.path.isfile(fmsvap.assignment_bam)
                assert mrp.rows == 1
        for mrp in split_merged_reads(mr, 2, included_reads='M'):
            fmsvaps = align_reads_to_ms_variations_part_as_list(mrp, padding, mss_version)
            for fmsvap in fmsvaps:
                assert fmsvap.merged_reads_part == mrp
                assert os.path.isfile(fmsvap.assignment_bam)
            assert mrp.rows == 2
    assert FullMSVMergedReadsPart.objects.filter(rows=1).count()
    assert FullMSVMergedReadsPart.objects.exclude(rows=1).count()
    assert FullMSVMergedReadsPart.objects.filter(rows=2).count()
    assert FullMSVMergedReadsPart.objects.exclude(rows=2).count()
    FullMSVMergedReadsPart.objects.all().delete()
    FullMSVariations.objects.all().delete()


@pytest.mark.django_db
def test_merge_aligned_merged_reads_parts(fmsv_merged_reads_d, fmsv_reads_fd, ac_134, requires_microsatellites):
    padding = 50
    mss_version = 0
    amplicon_collection = ac_134
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        amplicon_collection = mr.sample_reads.library.subclass.panel.amplicon_collection
        fmsv = get_full_ms_variations(amplicon_collection, padding, mss_version)
        fmsva_parts = []
        for mrp in split_merged_reads(mr, reads_chunk_size=1, included_reads='M'):
            fmsvaps = align_reads_to_ms_variations_part_as_list(mrp, padding, mss_version)
            for fmsvap in fmsvaps:
                assert fmsvap.ms_variations.id == fmsv.id
                assert fmsvap.merged_reads_part == mrp
                assert os.path.isfile(fmsvap.assignment_bam)
            fmsva_parts.append(fmsvaps)
        fmsvas = merge_fmsva_parts_as_list(fmsva_parts, reads_chunk_size=1, included_reads='M')
        for fmsva in fmsvas:
            assert fmsva.merged_reads == mr
            assert os.path.isfile(fmsva.sorted_assignment_bam)
            fmsva.delete()
        mr.fullmsvmergedreadspart_set.all().delete()
        fmsv.delete()


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
            fmsvas = list(align_reads_to_ms_variations(mr, padding, mss_version)) #included_reads=inc\
            for fmsva in fmsvas:
                for amp_id, histogram_reads in stream_group_alignemnts(fmsva):
                    generator_read_ids = []
                    for read_ids in histogram_reads.values():
                        generator_read_ids += read_ids
                    if inc == "M":
                        ref_reads_d = [seq.id for seq in fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp_id][RM]]
                    # TODO: inc == "F"
                    assert set(generator_read_ids) == set(ref_reads_d)
                fmsva.ms_variations.delete()
                fmsva.delete()
    msv.delete()

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db
def test_amplicons_mapping(fmsv_merged_reads_d, fmsv_reads_fd, full_ms_variations, amplicon_collection, requires_none_genotypes, write_her_file_flag):
    padding = 50
    mss_version = 0
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        sr = mr.sample_reads
        sr.write_her_files = write_her_file_flag
        sr.save()
        # for inc in ["M", "F"]:
        for inc in ["M"]:
            # test_ms_variations_index
            assert os.path.isfile(mr.assembled_fastq)
            fmsvas = list(align_reads_to_ms_variations(mr, padding, mss_version)) #included_reads=inc
            for fmsva in fmsvas:
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
                    if write_her_file_flag:
                        fmsvh_reads = {
                            'R1': itertools.chain(*[
                                SeqIO.parse(her.fastq1, "fastq") for her in
                                HistogramEntryReads.objects.filter(histogram=fmsvh)
                            ]),
                            'R2': itertools.chain(*[
                                SeqIO.parse(her.fastq2, "fastq") for her in
                                HistogramEntryReads.objects.filter(histogram=fmsvh)
                            ]),
                            'RM': itertools.chain(*[
                                SeqIO.parse(her.fastqm, "fastq") for her in
                                HistogramEntryReads.objects.filter(histogram=fmsvh)
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
                    if write_her_file_flag:
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
                #delete data
                fmsvv = fmsva.ms_variations
                fmsva.delete()
                assert not os.path.exists(fmsva.sorted_assignment_bam)
                fmsvv.delete()
                assert not os.path.exists(fmsvv.index_dump_dir)
            assert amps == set(fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED)) | \
                           (set(fmsv_reads_fd.keys(l_id, bc_id, UNASSEMBLED)) if \
                                inc == "F" else set())

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db
def test_genotype_mapping(fmsv_merged_reads_d, fmsv_reads_fd, full_ms_variations, amplicon_collection, requires_none_genotypes, write_her_file_flag):
    padding = 50
    mss_version = 0
    for (l_id, bc_id), mr in fmsv_merged_reads_d.items():
        sr = mr.sample_reads
        sr.write_her_files = write_her_file_flag
        sr.save()
        # for inc in ["M", "F"]:
        for inc in ["M"]:
            # test_ms_variations_index
            assert os.path.isfile(mr.assembled_fastq)
            fmsvas = list(align_reads_to_ms_variations(mr, padding, mss_version))  #included_reads=inc
            for fmsva in fmsvas:
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
                    if her.histogram.sample_reads.write_her_files:
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
                #delete data
                fmsvv = fmsva.ms_variations
                fmsva.delete()
                assert not os.path.exists(fmsva.sorted_assignment_bam)
                fmsvv.delete()
                assert not os.path.exists(fmsvv.index_dump_dir)
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


@pytest.mark.django_db
def test_separate_reads_by_genotypes_on_empty_fmsva(unknown_fmsv_merged_reads_d, unknown_fmsv_reads_fd, amplicon_collection, requires_none_genotypes):
    padding = 50
    mss_version = 0
    for (l_id, bc_id), mr in unknown_fmsv_merged_reads_d.items():
        # for inc in ["M", "F"]:
        for inc in ["M"]:
            # test_ms_variations_index
            assert os.path.isfile(mr.assembled_fastq)
            fmsvas = list(align_reads_to_ms_variations(mr, padding, mss_version))  #included_reads=inc
            for fmsva in fmsvas:
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
                        ref_reads_d = unknown_fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen]
                    else:  # inc == "F"
                        ref_reads_d = {
                            R1: unknown_fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen][R1] + \
                                unknown_fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp, gen][R1],
                            R2: unknown_fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen][R2] + \
                                unknown_fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp, gen][R2],
                            RM: unknown_fmsv_reads_fd[l_id, bc_id, ASSEMBLED, amp, gen][RM] + \
                                unknown_fmsv_reads_fd[l_id, bc_id, UNASSEMBLED, amp, gen][R1],
                        }
                    if her.histogram.sample_reads.write_her_files:
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
            assert gens == set(unknown_fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED, amp)) \
                           | set(unknown_fmsv_reads_fd.keys(l_id, bc_id, UNASSEMBLED, amp))
        else:
            ref_gens = set()
            for amp in unknown_fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED):
                for s in unknown_fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED, amp):
                    ref_gens.add(s)
            assert gens == ref_gens
            assert amps == set(unknown_fmsv_reads_fd.keys(l_id, bc_id, ASSEMBLED)) | \
                           (set(unknown_fmsv_reads_fd.keys(l_id, bc_id, UNASSEMBLED)) if \
                                inc == "F" else set())
        for fmsva in fmsvas:
            fmsvv = fmsva.ms_variations
            fmsva.delete()
            assert not os.path.exists(fmsva.sorted_assignment_bam)
            fmsvv.delete()
            assert not os.path.exists(fmsvv.index_dump_dir)

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db(transaction=True)
def test_run_mono_split_alignments(demultiplexing, sample_reads_d, fmsv_reads_fd, requires_amplicons, requires_microsatellites, requires_none_genotypes, write_her_file_flag):
    mss_version = 0
    ref_padding = 50
    reads_chunk_size = 1
    for sr in demultiplexing.samplereads_set.all():
        sr.write_her_files = write_her_file_flag
        sr.save()
        amplicon_collection = sr.library.subclass.panel.amplicon_collection
        msv = get_full_ms_variations(amplicon_collection, ref_padding, mss_version)
    herss = {inc: set() for inc in ["M"]}  # TODO: "F"
    for inc in herss.keys():
        # FIXME?
        merged_reads = [merge(sr) for sr in demultiplexing.samplereads_set.all()]
        fmsv_merged_reads_parts_lists = [split_merged_reads_as_list(mr, reads_chunk_size, inc) for mr in merged_reads]
        fmsva_parts_lists = [list(map(align_reads_to_ms_variations_part_as_list, fmsv_merged_reads_parts_list, itertools.repeat(ref_padding), itertools.repeat(mss_version))) for fmsv_merged_reads_parts_list in fmsv_merged_reads_parts_lists]
        fmsvas = [merge_fmsva_parts(fmsva_parts_list, reads_chunk_size, inc) for fmsva_parts_list in fmsva_parts_lists]
        fhers_list = [list(separate_reads_by_genotypes(fmsva)) for fmsva in fmsvas]
        for fhers in fhers_list:
            hers_gen = fhers
            for her in hers_gen:
                herss[inc].add(her)
    for inc, hers in herss.items():
        parts = set()
        for her in hers:
            amp = her.histogram.amplicon_id
            sr = her.histogram.sample_reads
            bc = sr.barcoded_content_id
            l_id = sr.library_id
            gen = frozenset((msg.microsatellite_id, msg.repeat_number) for \
                msg in her.microsatellite_genotypes.genotypes)
            parts.add((l_id, bc, amp, gen))

            her_fnames_d = {
                R1: her.fastq1,
                R2: her.fastq2,
                RM: her.fastqm,
            }

            if inc == "M":
                ref_reads_d = fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen]
            else:  # inc == "F"
                ref_reads_d = {
                    R1: fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen][R1] + \
                        fmsv_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R1],
                    R2: fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen][R2] + \
                        fmsv_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R2],
                    RM: fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen][RM] + \
                        fmsv_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R1],
                }
            if her.histogram.sample_reads.write_her_files:
                for r in [R1, R2, RM]:
                    assert set(srs_to_tups(  # TODO: get informative error on genotyping mismatch
                        SeqIO.parse(her_fnames_d[r], "fastq"))
                    ) == \
                    set(srs_to_tups(
                        ref_reads_d[r]
                    ))
            assert her.num_reads == \
                len(ref_reads_d[RM])
        ref_parts = set()
        for key_tups in fmsv_reads_fd:
            if len(key_tups) < 5:
                continue
            l_id2, bc2, t2, amp2, gen2 = key_tups
            if t2 == ASSEMBLED or inc == "F":
                ref_parts.add((l_id2, bc2, amp2, gen2))
        assert ref_parts == parts

    for Model in [
        FullMSVHistogram,
        FullMSVariations,  # FIXME: get these from outside.
        FullMSVAssignment,
        FullMSVAssignmentPart,
        FullMSVMergedReads,
        FullMSVMergedReadsPart,
        HistogramEntryReads,  # FIXME: get these from outside.
    ]:
        Model.objects.all().delete()

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db(transaction=True)
def test_run_mono_split_alignments_split_mapping(demultiplexing, sample_reads_d, fmsv_reads_fd, requires_amplicons, requires_microsatellites, requires_none_genotypes, write_her_file_flag):
    mss_version = 0
    ref_padding = 50
    reads_chunk_size = 1
    amplicons_chunk_size = 1
    for sr in demultiplexing.samplereads_set.all():
        sr.write_her_files = write_her_file_flag
        sr.save()
        total_amplicon_collection = sr.library.subclass.panel.amplicon_collection
        all_amplicons = total_amplicon_collection.amplicons.order_by('id')
        amplicons_splitted = grouper(amplicons_chunk_size,
                                     all_amplicons)  # split the amplicons to create smaller amplicon collections
        for amplicon_subgroup in amplicons_splitted:
            amplicon_collection = AmpliconCollection.custom_get_or_create(amplicons=amplicon_subgroup)
            msv = get_full_ms_variations(amplicon_collection, ref_padding, mss_version)
    herss = {inc: set() for inc in ["M"]}  # TODO: "F"
    for inc in herss.keys():
        # FIXME?
        merged_reads = [merge(sr) for sr in demultiplexing.samplereads_set.all()]
        fmsv_merged_reads_parts_lists = [split_merged_reads_as_list(mr, reads_chunk_size, inc) for mr in merged_reads]
        fmsva_parts_lists = [
            list(
                map(
                    align_reads_to_ms_variations_part_as_list,
                        fmsv_merged_reads_parts_list,
                        itertools.repeat(ref_padding),
                        itertools.repeat(mss_version),
                        itertools.repeat(amplicons_chunk_size),
                )
            ) for fmsv_merged_reads_parts_list in fmsv_merged_reads_parts_lists
        ]
        fmsvass = [list(merge_fmsva_parts(fmsva_parts_list, reads_chunk_size, inc)) for fmsva_parts_list in fmsva_parts_lists]
        fhers_list = [list(separate_reads_by_genotypes(fmsva)) for fmsvas in fmsvass for fmsva in fmsvas]
        for fhers in fhers_list:
            hers_gen = fhers
            for her in hers_gen:
                herss[inc].add(her)
    for inc, hers in herss.items():
        parts = set()
        for her in hers:
            amp = her.histogram.amplicon_id
            sr = her.histogram.sample_reads
            bc = sr.barcoded_content_id
            l_id = sr.library_id
            gen = frozenset((msg.microsatellite_id, msg.repeat_number) for \
                msg in her.microsatellite_genotypes.genotypes)
            parts.add((l_id, bc, amp, gen))

            her_fnames_d = {
                R1: her.fastq1,
                R2: her.fastq2,
                RM: her.fastqm,
            }

            if inc == "M":
                ref_reads_d = fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen]
            else:  # inc == "F"
                ref_reads_d = {
                    R1: fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen][R1] + \
                        fmsv_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R1],
                    R2: fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen][R2] + \
                        fmsv_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R2],
                    RM: fmsv_reads_fd[l_id, bc, ASSEMBLED, amp, gen][RM] + \
                        fmsv_reads_fd[l_id, bc, UNASSEMBLED, amp, gen][R1],
                }
            if her.histogram.sample_reads.write_her_files:
                for r in [R1, R2, RM]:
                    assert set(srs_to_tups(  # TODO: get informative error on genotyping mismatch
                        SeqIO.parse(her_fnames_d[r], "fastq"))
                    ) == \
                    set(srs_to_tups(
                        ref_reads_d[r]
                    ))
            assert her.num_reads == \
                len(ref_reads_d[RM])
        ref_parts = set()
        for key_tups in fmsv_reads_fd:
            if len(key_tups) < 5:
                continue
            l_id2, bc2, t2, amp2, gen2 = key_tups
            if t2 == ASSEMBLED or inc == "F":
                ref_parts.add((l_id2, bc2, amp2, gen2))
        assert ref_parts == parts

    for Model in [
        FullMSVHistogram,
        FullMSVariations,  # FIXME: get these from outside.
        FullMSVAssignment,
        FullMSVAssignmentPart,
        FullMSVMergedReads,
        FullMSVMergedReadsPart,
        HistogramEntryReads,  # FIXME: get these from outside.
    ]:
        Model.objects.all().delete()
