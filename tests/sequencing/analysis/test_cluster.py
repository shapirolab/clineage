import pytest
import os
from Bio import SeqIO

from sequencing.analysis.adamiya import merge

from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM, \
    srs_to_tups, rc_srs_to_tups
from tests.sequencing.analysis.reads_dict import ASSEMBLED, UNASSEMBLED

from distributed.utils_test import gen_cluster
from distributed.executor import Future, as_completed


@pytest.mark.django_db(transaction=True)
def test_remote_runmerge_single(adam_reads_fd, sample_reads_d):
    @gen_cluster(ncores=[('127.0.0.1', 1)], executor=True)
    def _inner_test_remote_runmerge(e, s, *args):
        for (l_id, bc_id), sr in sample_reads_d.items():
            mr_f = e.submit(merge, sr)
            assert isinstance(mr_f, Future)
            mr = yield mr_f._result()
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
    _inner_test_remote_runmerge()


@pytest.mark.django_db(transaction=True)
def test_remote_runmerge_parallel(adam_reads_fd, sample_reads_d):
    @gen_cluster(ncores=[('127.0.0.1', 1)]*4, executor=True)
    def _inner_test_remote_runmerge(e, s, *args):
        futures = []
        for (l_id, bc_id), sr in sample_reads_d.items():
            mr_f = e.submit(merge, sr)
            futures.append((l_id, bc_id, sr, mr_f))
        for (l_id, bc_id, sr, mr_f) in futures:
            assert isinstance(mr_f, Future)
            mr = yield mr_f._result()
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
    _inner_test_remote_runmerge()


@pytest.mark.django_db(transaction=True)
def test_map_runmerge(adam_reads_fd, sample_reads_d):
    @gen_cluster(executor=True)
    def _inner_test_map_runmerge(e, s, a, b):
        def _inner_merge(key_tup, sr):
            mr = merge(sr)
            return key_tup, mr
        fs = e.map(_inner_merge, *zip(*list(sample_reads_d.items())))
        for f in as_completed(iter(fs)):
            key_tup, mr = yield f._result()
            l_id, bc_id = key_tup
            assert set(srs_to_tups(SeqIO.parse(mr.assembled_fastq, "fastq"))) == \
                set(srs_to_tups(adam_reads_fd[l_id, bc_id, ASSEMBLED][RM]))
            assert set(srs_to_tups(SeqIO.parse(mr.discarded_fastq, "fastq"))) == set()
            assert set(srs_to_tups(SeqIO.parse(mr.unassembled_forward_fastq, "fastq"))) == \
                set(srs_to_tups(adam_reads_fd[l_id, bc_id, UNASSEMBLED][R1]))
            assert set(srs_to_tups(SeqIO.parse(mr.unassembled_reverse_fastq, "fastq"))) == \
                set(rc_srs_to_tups(adam_reads_fd[l_id, bc_id, UNASSEMBLED][R2]))
            mr.delete()
    _inner_test_map_runmerge()
