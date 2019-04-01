import pytest
import os
from Bio import SeqIO
from frogress import bar
from sequencing.analysis.full_msv.full_msv import merge
from sequencing.analysis.full_msv.parallel_fmsv import run_parallel, run_parallel_split_alignments

from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM, \
    srs_to_tups, rc_srs_to_tups
from tests.sequencing.analysis.full_msv.reads_dict import ASSEMBLED, UNASSEMBLED


from distributed import as_completed
from distributed.client import Future
from sequencing.analysis.full_msv.models import FullMSVMergedReads,\
    FullMSVariations, FullMSVHistogram, FullMSVAssignment
from sequencing.analysis.models import HistogramEntryReads, MicrosatelliteHistogramGenotype
from sequencing.analysis.full_msv.full_msv import get_full_ms_variations


@pytest.mark.django_db(transaction=True)
def test_remote_runmerge_single(executor, fmsv_reads_fd, sample_reads_d):
    for (l_id, bc_id), sr in sample_reads_d.items():
        mr_f = executor.submit(merge, sr)
        assert isinstance(mr_f, Future)
        mr = mr_f.result()
        assert mr.sample_reads.id == sr.id
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


@pytest.mark.django_db(transaction=True)
def test_remote_runmerge_parallel(executor, fmsv_reads_fd, sample_reads_d):
    futures = []
    for (l_id, bc_id), sr in sample_reads_d.items():
        mr_f = executor.submit(merge, sr)
        futures.append((l_id, bc_id, sr, mr_f))
    for (l_id, bc_id, sr, mr_f) in futures:
        assert isinstance(mr_f, Future)
        mr = mr_f.result()
        assert mr.sample_reads.id == sr.id
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


@pytest.mark.django_db(transaction=True)
def test_map_runmerge(executor, fmsv_reads_fd, sample_reads_d):
    def _inner_merge(key_tup, sr):
        mr = merge(sr)
        return key_tup, mr
    fs = executor.map(_inner_merge, *zip(*list(sample_reads_d.items())))
    for f in as_completed(iter(fs)):
        key_tup, mr = f.result()
        l_id, bc_id = key_tup
        assert set(srs_to_tups(SeqIO.parse(mr.assembled_fastq, "fastq"))) == \
            set(srs_to_tups(fmsv_reads_fd[l_id, bc_id, ASSEMBLED][RM]))
        assert set(srs_to_tups(SeqIO.parse(mr.discarded_fastq, "fastq"))) == set()
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_forward_fastq, "fastq"))) == \
            set(srs_to_tups(fmsv_reads_fd[l_id, bc_id, UNASSEMBLED][R1]))
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_reverse_fastq, "fastq"))) == \
            set(rc_srs_to_tups(fmsv_reads_fd[l_id, bc_id, UNASSEMBLED][R2]))
        mr.delete()

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db(transaction=True)
def test_run_parallel(executor, demultiplexing, sample_reads_d, fmsv_reads_fd, requires_amplicons, requires_microsatellites, requires_none_genotypes, write_her_file_flag):
    mss_version = 0
    ref_padding = 50
    for sr in demultiplexing.samplereads_set.all():
        sr.write_her_files = write_her_file_flag
        sr.save()
        amplicon_collection = sr.library.subclass.panel.amplicon_collection
        msv = get_full_ms_variations(amplicon_collection, ref_padding, mss_version)
    herss = {inc: set() for inc in ["M"]}  # TODO: "F"
    for inc in herss.keys():
        # FIXME?
        res = run_parallel(executor, demultiplexing.samplereads_set.all(), inc)
        merged_reads = next(res)
        fmsvas = next(res)
        fhers_list = next(res)
        for fhers in as_completed(fhers_list):
            hers_gen = fhers.result()
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
        FullMSVMergedReads,
        HistogramEntryReads,  # FIXME: get these from outside.
    ]:
        Model.objects.all().delete()

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db(transaction=True)
def test_run_parallel_small_size_amplicon_collection(executor, demultiplexing, sample_reads_d, fmsv_reads_fd, requires_amplicons, requires_microsatellites, requires_none_genotypes, write_her_file_flag):
    mss_version = 0
    ref_padding = 50
    for sr in demultiplexing.samplereads_set.all():
        sr.write_her_files = write_her_file_flag
        sr.save()
        amplicon_collection = sr.library.subclass.panel.amplicon_collection
        msv = get_full_ms_variations(amplicon_collection, ref_padding, mss_version)
    herss = {inc: set() for inc in ["M"]}  # TODO: "F"
    for inc in herss.keys():
        # FIXME?
        res = run_parallel(executor, demultiplexing.samplereads_set.all(), inc, amp_col_size=1)
        merged_reads = next(res)
        fmsvas = next(res)
        fhers_list = next(res)
        for fhers in as_completed(fhers_list):
            hers_gen = fhers.result()
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
        FullMSVMergedReads,
        HistogramEntryReads,  # FIXME: get these from outside.
    ]:
        Model.objects.all().delete()

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db(transaction=True)
def test_run_parallel_small_size_amplicon_collection(executor, demultiplexing, sample_reads_d, fmsv_reads_fd, requires_amplicons, requires_microsatellites, requires_none_genotypes, write_her_file_flag):
    mss_version = 0
    ref_padding = 50
    for sr in demultiplexing.samplereads_set.all():
        sr.write_her_files = write_her_file_flag
        sr.save()
        amplicon_collection = sr.library.subclass.panel.amplicon_collection
        msv = get_full_ms_variations(amplicon_collection, ref_padding, mss_version)
    herss = {inc: set() for inc in ["M"]}  # TODO: "F"
    for inc in herss.keys():
        # FIXME?
        res = run_parallel(executor, demultiplexing.samplereads_set.all(), inc, amp_col_size=1)
        merged_reads = next(res)
        fmsvas = next(res)
        fhers_list = next(res)
        for fhers in as_completed(fhers_list):
            hers_gen = fhers.result()
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
        FullMSVMergedReads,
        HistogramEntryReads,  # FIXME: get these from outside.
    ]:
        Model.objects.all().delete()

@pytest.mark.parametrize('write_her_file_flag', [True, False])
@pytest.mark.django_db(transaction=True)
def test_run_parallel_split_alignments(executor, demultiplexing, sample_reads_d, fmsv_reads_fd, requires_amplicons, requires_microsatellites, requires_none_genotypes, write_her_file_flag):
    mss_version = 0
    ref_padding = 50
    for sr in demultiplexing.samplereads_set.all():
        sr.write_her_files = write_her_file_flag
        sr.save()
        amplicon_collection = sr.library.subclass.panel.amplicon_collection
        msv = get_full_ms_variations(amplicon_collection, ref_padding, mss_version)
    herss = {inc: set() for inc in ["M"]}  # TODO: "F"
    for inc in herss.keys():
        # FIXME?
        rpsa_gen = run_parallel_split_alignments(executor, demultiplexing.samplereads_set.all(), inc)
        merged_reads = next(rpsa_gen)
        fmsv_merged_reads_parts_lists = next(rpsa_gen)
        fmsva_parts_lists = next(rpsa_gen)
        merged_fmsvas = next(rpsa_gen)
        fhers_list = next(rpsa_gen)
        for fhers in fhers_list:
            hers_gen = fhers.result()
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
        FullMSVMergedReads,
        HistogramEntryReads,  # FIXME: get these from outside.
    ]:
        Model.objects.all().delete()
