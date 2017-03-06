import pytest
import os
from Bio import SeqIO
import vcf

from sequencing.analysis.full_msv.models import SampleReads
from sequencing.analysis.models_common import SNPHistogramGenotype, MicrosatelliteHistogramGenotype

from sequencing.analysis.snps.models import ReadsAlignment, VCFReads
from sequencing.analysis.snps.snp_calling import sam_to_bam, sort_and_index_bam_file, get_alignment_reference, \
    make_faidx, mpileup
from tests.sequencing.runs.conftest import *
from tests.sequencing.analysis.snps.amp_fasta import VARS_AMP_1234, SAM_EXAMPLE_TXT
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.amplicons.conftest import *
from tests.targeted_enrichment.planning.conftest import *
from tests.sequencing.analysis.reads_dict_tools import R1, R2, RM, NUM_READS
from tests.flat_dict import FlatDict
from misc.utils import unique_file_cm
from tests.sequencing.analysis.snps.reads_dict import READS_DICT_FULL_SNP


@pytest.fixture()
def requires_none_genotypes(request, transactional_db):
    MicrosatelliteHistogramGenotype.objects.get_or_create(microsatellite=None,
        defaults=dict(repeat_number=1))
    SNPHistogramGenotype.objects.get_or_create(snp=None, defaults=dict(base=""))


@pytest.fixture(scope="session")
def snp_reads_fd():
    return FlatDict(READS_DICT_FULL_SNP, [R1, R2, RM])


@pytest.fixture()
def amp_1234_full_msv():
    return VARS_AMP_1234


@pytest.fixture()
def sam_example():
    sam_file = get_unique_path("sam")
    with open(sam_file, 'w') as f:
        for text in SAM_EXAMPLE_TXT:
            for line in text:
                f.write(line+'\n')
    return sam_file


@pytest.yield_fixture(scope="session")
def sample_reads_files_d(snp_reads_fd, temp_storage):
    d = {}
    for l_id, l_d in snp_reads_fd.items():
        for bc_id, r_d in l_d.reads():
            fastq_r1 = get_unique_path("fastq")
            fastq_r2 = get_unique_path("fastq")
            SeqIO.write(r_d[R1], fastq_r1, "fastq")
            SeqIO.write(r_d[R2], fastq_r2, "fastq")
            d[l_id, bc_id] = {R1: fastq_r1, R2: fastq_r2, NUM_READS: len(r_d[R1])}
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
            num_reads=f_d[NUM_READS],
        )
        # So our objects don't have "special" objects in fields
        sr = SampleReads.objects.get(pk=sr.pk)
        d[l_id, bc_id] = sr
    yield d
    for sr in d.values():
        os.unlink(sr.fastq1)
        os.unlink(sr.fastq2)


@pytest.yield_fixture()
def snp_reads_d(sample_reads_d, snp_reads_fd, temp_storage):
    d = {}
    for (l_id, bc_id), sr in sample_reads_d.items():
        d[l_id, bc_id] = snp_reads_fd[(l_id, bc_id)]
    yield d



@pytest.yield_fixture()
def alignment_ref(snp_1234, temp_storage):
    acbwai = get_alignment_reference(snp_1234)
    yield acbwai

    acbwai.delete()
    # for file in acbwai.files:
    #     os.unlink(file)


@pytest.yield_fixture()
def read_alignment(sample_reads_d, alignment_ref, sam_example, temp_storage):
    d = {}
    bam_file = sam_to_bam(sam_example)
    bam_sort_file = sort_and_index_bam_file(bam_file)

    for (l_id, bc_id), sr in sample_reads_d.items():
        ra = ReadsAlignment.objects.create(
            sample_read=sr,
            alignment_reference=alignment_ref,
            bam_file=bam_sort_file,
        )
        # So our objects don't have "special" objects in fields
        ra = ReadsAlignment.objects.get(pk=ra.pk)
        d[l_id, bc_id] = ra
    yield d

    for ra in d.values():
        ra.delete()


@pytest.yield_fixture()
def vcf_object(read_alignment, requires_snp_targets, requires_snps, temp_storage):
    d = {}
    for (l_id, bc_id), ra in read_alignment.items():
        faidx_file = make_faidx(ra.alignment_reference.fasta_file)
        vcf_file = mpileup(ra.bam_file, ra.alignment_reference.fasta_file)
        vcf = VCFReads.objects.create(
            reads_alignment=ra,
            vcf_file=vcf_file,
        )

        os.unlink(faidx_file)
        assert not os.path.exists(faidx_file)

        # So our objects don't have "special" objects in fields
        vcf = VCFReads.objects.get(pk=vcf.pk)
        d[l_id, bc_id] = vcf
    yield d

    for vcf in d.values():
        os.unlink(vcf.vcf_file)


@pytest.fixture()
def vcf_rec(vcf_object):
    for (l_id, bc_id), vc in vcf_object.items():
        vcf_reader = vcf.Reader(open(vc.vcf_file, 'r'))
        for row in vcf_reader:
            if row.CHROM == '7':
                return row
