import pytest
import os
from Bio import SeqIO

from sequencing.analysis.adamiya import merge, run_parallel

from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM, \
    srs_to_tups, rc_srs_to_tups
from tests.sequencing.analysis.reads_dict import ASSEMBLED, UNASSEMBLED

from distributed.utils_test import cluster
from distributed.executor import Future, as_completed, Executor

from sequencing.analysis.models import AdamMergedReads, AdamReadsIndex, \
    AdamMarginAssignment, AdamAmpliconReads, AdamHistogram, AdamMSVariations, \
    HistogramEntryReads, MicrosatelliteHistogramGenotype

import django

@pytest.yield_fixture(scope="session")
def executor(): 
    with cluster(4) as (d, workers):
        yield Executor(("127.0.0.1", d["port"]))


@pytest.mark.django_db(transaction=True)
def test_remote_runmerge_single(executor, adam_reads_fd, sample_reads_d):
    for (l_id, bc_id), sr in sample_reads_d.items():
        mr_f = executor.submit(merge, sr)
        assert isinstance(mr_f, Future)
        mr = mr_f.result()
        assert mr.sample_reads.id == sr.id
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


@pytest.mark.django_db(transaction=True)
def test_remote_runmerge_parallel(executor, adam_reads_fd, sample_reads_d):
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


@pytest.mark.django_db(transaction=True)
def test_map_runmerge(executor, adam_reads_fd, sample_reads_d):
    def _inner_merge(key_tup, sr):
        mr = merge(sr)
        return key_tup, mr
    fs = executor.map(_inner_merge, *zip(*list(sample_reads_d.items())))
    for f in as_completed(iter(fs)):
        key_tup, mr = f.result()
        l_id, bc_id = key_tup
        assert set(srs_to_tups(SeqIO.parse(mr.assembled_fastq, "fastq"))) == \
            set(srs_to_tups(adam_reads_fd[l_id, bc_id, ASSEMBLED][RM]))
        assert set(srs_to_tups(SeqIO.parse(mr.discarded_fastq, "fastq"))) == set()
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_forward_fastq, "fastq"))) == \
            set(srs_to_tups(adam_reads_fd[l_id, bc_id, UNASSEMBLED][R1]))
        assert set(srs_to_tups(SeqIO.parse(mr.unassembled_reverse_fastq, "fastq"))) == \
            set(rc_srs_to_tups(adam_reads_fd[l_id, bc_id, UNASSEMBLED][R2]))
        mr.delete()


@pytest.mark.xfail(django.db.connection.vendor == "sqlite",
    raises=MicrosatelliteHistogramGenotype.MultipleObjectsReturned,
    reason="SQLite isn't transactional enough for distributing.")
@pytest.mark.django_db(transaction=True)
def test_run_parallel(executor, demultiplexing, sample_reads_d, adam_reads_fd, requires_amplicons, requires_microsatellites):
    herss = {inc: set(run_parallel(executor, demultiplexing, inc)) for \
        inc in ["M", "F"]}
    for inc, hers in herss.items():
        parts = set()
        for her in hers:
            amp = her.amplicon_id
            sr = her.histogram.subclass.amplicon_reads.margin_assignment.reads_index.merged_reads.sample_reads
            bc = sr.barcoded_content_id
            l_id = sr.library_id
            gen = frozenset((msg.microsatellite_id, msg.repeat_number) for \
                msg in her.microsatellite_genotypes.all())
            parts.add((l_id, bc, amp, gen))
            
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
        ref_parts = set()
        for key_tups in adam_reads_fd:
            if len(key_tups) < 5:
                continue
            l_id2, bc2, t2, amp2, gen2 = key_tups
            if t2 == ASSEMBLED or inc == "F":
                ref_parts.add((l_id2, bc2, amp2, gen2))
        assert ref_parts == parts
    for Model in [
        AdamHistogram,
        AdamMSVariations,  # FIXME: get these from outside.
        AdamAmpliconReads,
        AdamMarginAssignment,
        AdamReadsIndex,
        AdamMergedReads,
    ]:
        Model.objects.all().delete()
