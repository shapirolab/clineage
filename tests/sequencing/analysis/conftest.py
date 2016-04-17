import pytest
import os
import uuid
from django.conf import settings
import itertools
from Bio import SeqIO

from sequencing.analysis.models import SampleReads, AdamMergedReads, \
    AdamReadsIndex, AdamMarginAssignment, AdamAmpliconReads, \
    AdamMSVariations, AdamHistogram
from targeted_enrichment.amplicons.models import Amplicon
from sequencing.analysis.models import LEFT, RIGHT
from misc.utils import get_unique_path

from tests.sequencing.runs.conftest import *
from tests.sequencing.analysis.pu_28727_adam_ms_variations import VARS_28727
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.amplicons.conftest import *
from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM
from tests.flat_dict import FlatDict
from tests.sequencing.analysis.reads_dict import READS_DICT_ADAM, ASSEMBLED, \
    UNASSEMBLED


file_fixtures_path = os.path.join(*(os.path.split(os.path.dirname(os.path.realpath(__file__)))[:-1] + ("ngs_fixtures",)))


@pytest.fixture(scope="session")
def adam_reads_fd():
    return FlatDict(READS_DICT_ADAM, [R1, R2, RM])


@pytest.yield_fixture(scope="session")
def sample_reads_files_d(adam_reads_fd, temp_storage):
    d = {}
    for l_id, l_d in adam_reads_fd.items():
        for bc_id, r_d in l_d.reads():
            fastq_r1 = get_unique_path("fastq")
            fastq_r2 = get_unique_path("fastq")
            SeqIO.write(r_d[R1], fastq_r1, "fastq")
            SeqIO.write(r_d[R2], fastq_r2, "fastq")
            d[l_id, bc_id] = {R1: fastq_r1, R2: fastq_r2}
    yield d
    for f_d in d.values():
        os.unlink(f_d[R1])
        os.unlink(f_d[R2])


@pytest.yield_fixture()
def sample_reads_d(sample_reads_files_d, demultiplexing, require_magicals):
    d = {}
    for (l_id, bc_id), f_d in sample_reads_files_d.items():
        fastq_r1 = get_unique_path("fastq")
        fastq_r2 = get_unique_path("fastq")
        os.symlink(f_d[R1], fastq_r1)
        os.symlink(f_d[R2], fastq_r2)
        sr = SampleReads.objects.create(
            demux=demultiplexing,
            barcoded_content=MagicalPCR1BarcodedContent.objects.get(id=bc_id),
            library=MagicalPCR1Library.objects.get(id=l_id),
            fastq1=fastq_r1,
            fastq2=fastq_r2,
        )
        d[l_id, bc_id] = sr
    yield d
    for sr in d.values():
        # sr.delete()
        os.unlink(sr.fastq1)
        os.unlink(sr.fastq2)


@pytest.yield_fixture(scope="session")
def adam_merged_reads_files_d(adam_reads_fd, temp_storage):
    d = {}
    for l_id, l_d in adam_reads_fd.items():
        for bc_id, s_d in l_d.items():
            assembled_fastq = get_unique_path("fastq")
            unassembled_forward_fastq = get_unique_path("fastq")
            unassembled_reverse_fastq = get_unique_path("fastq")
            discarded_fastq = get_unique_path("fastq")
            SeqIO.write(s_d[ASSEMBLED][RM], assembled_fastq, "fastq")
            SeqIO.write(s_d[UNASSEMBLED][R1], unassembled_forward_fastq, "fastq")
            SeqIO.write(s_d[UNASSEMBLED][R2], unassembled_reverse_fastq, "fastq")
            SeqIO.write((), discarded_fastq, "fastq")
            d[l_id, bc_id] = {
                ASSEMBLED: assembled_fastq,
                (UNASSEMBLED, R1): unassembled_forward_fastq,
                (UNASSEMBLED, R2): unassembled_reverse_fastq,
                None: discarded_fastq
            }
    yield d
    for f_d in d.values():
        os.unlink(f_d[ASSEMBLED])
        os.unlink(f_d[UNASSEMBLED, R1])
        os.unlink(f_d[UNASSEMBLED, R2])
        os.unlink(f_d[None])


@pytest.yield_fixture()
def adam_merged_reads_d(adam_merged_reads_files_d, sample_reads_d):
    d = {}
    for (l_id, bc_id), f_d in adam_merged_reads_files_d.items():
        dst_prefix = get_unique_path()
        os.symlink(
            f_d[ASSEMBLED],
            "{}.assembled.fastq".format(dst_prefix)
        )
        os.symlink(
            f_d[None],
            "{}.discarded.fastq".format(dst_prefix)
        )
        os.symlink(
            f_d[UNASSEMBLED, R1],
            "{}.unassembled.forward.fastq".format(dst_prefix)
        )
        os.symlink(
            f_d[UNASSEMBLED, R2],
            "{}.unassembled.reverse.fastq".format(dst_prefix)
        )
        mr = AdamMergedReads.objects.create(
            sample_reads=sample_reads_d[l_id, bc_id],
            assembled_fastq="{}.assembled.fastq".format(dst_prefix),
            discarded_fastq="{}.discarded.fastq".format(dst_prefix),
            unassembled_forward_fastq="{}.unassembled.forward.fastq".format(dst_prefix),
            unassembled_reverse_fastq="{}.unassembled.reverse.fastq".format(dst_prefix),
        )
        d[l_id, bc_id] = mr
    yield d
    for mr in d.values():
        # mr.delete()
        os.unlink(mr.assembled_fastq)
        os.unlink(mr.discarded_fastq)
        os.unlink(mr.unassembled_forward_fastq)
        os.unlink(mr.unassembled_reverse_fastq)


@pytest.yield_fixture()
def adam_amplicon_reads_files_d(adam_reads_fd, temp_storage):
    d = {}
    for l_id, l_d in adam_reads_fd.items():
        for bc_id, s_d in l_d.items():
            # M
            d[l_id, bc_id, "M"] = {}
            for amp, r_d in s_d.sub(ASSEMBLED).reads():
                f_d = {
                    RM: get_unique_path("fastq"),
                    R1: get_unique_path("fastq"),
                    R2: get_unique_path("fastq"),
                }
                for r in [R1, R2, RM]:
                    SeqIO.write(r_d[r], f_d[r], "fastq")
                d[l_id, bc_id, "M"][amp] = f_d
            # F
            d[l_id, bc_id, "F"] = {}
            all_amps = set(s_d.sub(ASSEMBLED).keys()) | set(s_d.sub(UNASSEMBLED).keys())
            for amp in all_amps:
                # Assuming there are no amplicons only present in UNASSEMBLED
                f_d = {
                    RM: get_unique_path("fastq"),
                    R1: get_unique_path("fastq"),
                    R2: get_unique_path("fastq"),
                }
                for r in [R1, R2]:
                    SeqIO.write(itertools.chain(
                        s_d[ASSEMBLED, amp][r],
                        s_d[UNASSEMBLED, amp][r],
                        ), f_d[r], "fastq")
                SeqIO.write(itertools.chain(
                    s_d[ASSEMBLED, amp][RM],
                    s_d[UNASSEMBLED, amp][R1],
                    ), f_d[RM], "fastq")
                d[l_id, bc_id, "F"][amp] = f_d
    yield d
    for f_d_d in d.values():
        for f_d in f_d_d.values():
            os.unlink(f_d[RM])
            os.unlink(f_d[R1])
            os.unlink(f_d[R2])


@pytest.yield_fixture()
def _chain(adam_amplicon_reads_files_d, adam_merged_reads_d):
    d = {}
    for (l_id, bc_id, inc), f_d_d in adam_amplicon_reads_files_d.items():
        dst_dir = get_unique_path()
        os.mkdir(dst_dir)
        ari = AdamReadsIndex.objects.create(
            merged_reads=adam_merged_reads_d[l_id, bc_id],
            included_reads=inc,  # Merged and unassembled_forward
            index_dump_dir=dst_dir,
            padding=5,
        )
        fake_sam = get_unique_path("sam")
        with open(fake_sam, "wb") as f:
            pass
        ama = AdamMarginAssignment.objects.create(
            reads_index=ari,
            assignment_sam=fake_sam,
        )
        d[l_id, bc_id, inc] = ari, ama
    yield d
    for ari, ama in d.values():
        os.rmdir(ari.index_dump_dir)
        os.unlink(ama.assignment_sam)


@pytest.yield_fixture()
def adam_amplicon_reads_d(adam_amplicon_reads_files_d, _chain, requires_amplicons):
    d = {}
    extra_dirs = []
    extra_files = []
    for (l_id, bc_id, inc), f_d_d in adam_amplicon_reads_files_d.items():
        ari, ama = _chain[l_id, bc_id, inc]
        for amp, f_d in f_d_d.items():
            f_d2 = {}
            for r in [R1, R2, RM]:
                f_d2[r] = get_unique_path("fastq")
                os.symlink(f_d[r], f_d2[r])
            aar = AdamAmpliconReads.objects.create(
                margin_assignment=ama,
                amplicon_id=amp,
                fastqm=f_d2[RM],
                fastq1=f_d2[R1],
                fastq2=f_d2[R2],
            )
            d[l_id, bc_id, inc, amp] = aar
    yield d
    for aar in d.values():
        # aar.margin_assignment.reads_index.delete()
        # aar.margin_assignment.delete()
        # aar.delete()
        os.unlink(aar.fastqm)
        os.unlink(aar.fastq1)
        os.unlink(aar.fastq2)
