import pytest
import os
import io
from Bio import SeqIO
import itertools

from sequencing.analysis.snps.snp_calling import get_fasta_file, get_alignment_reference, align_reads, pileup, \
    check_targets_id, parse_vcf
from targeted_enrichment.planning.models import SNP
from sequencing.analysis.snps.parse_snps import retrieve_explicit_snps_positions
from sequencing.analysis.snps.models import ReadsAlignment, VCFReads, SNPReads, _parse_vcf_file
from sequencing.analysis.models_common import SampleReads
from tests.sequencing.analysis.reads_dict_tools import strip_fasta_records


@pytest.mark.django_db
def test_create_fasta_file(amp_1234_full_msv, snp_1234):
    with get_fasta_file(snp_1234) as fasta:
        assert os.path.exists(fasta)
        variations = set(strip_fasta_records(SeqIO.parse(fasta, "fasta")))
        os.unlink(fasta)
    assert variations == amp_1234_full_msv



@pytest.mark.django_db
def test_build_bwa_index(amp_1234_full_msv, snp_1234):

    fmsv = get_alignment_reference(snp_1234)
    fasta = fmsv.fasta_file
    variations = set(strip_fasta_records(SeqIO.parse(fasta, "fasta")))
    assert variations == amp_1234_full_msv
    os.unlink(fmsv.fasta_file)
    fmsv.delete()
    for fp in fmsv.files:
        assert not os.path.exists(fp)
    assert not os.path.exists(fmsv.fasta_file)
    assert not os.path.exists(fmsv.index_dump_dir)


@pytest.mark.django_db
def test_snp_alignment(sample_reads_d, snp_reads_fd, snp_1234):
    assert not sample_reads_d.items() == []
    # acbwai = get_alignment_reference(snp_1234)
    srs = set()
    snp_srs = set()
    for (l_id, bc_id), sr in sample_reads_d.items():
        # test_snp_alignment
        acbwai = get_alignment_reference(snp_1234)  # AmpliconCollectionBWAIndex object
        reads_alignment = align_reads(sr, acbwai)  # ReadsAlignment object
        assert reads_alignment.sample_read == sr

        # test_reads for snps
        assert ReadsAlignment.objects.filter(
            sample_read=reads_alignment.sample_read,
        ).count() > 0
        reads_alignment_validation = ReadsAlignment.objects.get(
            sample_read=reads_alignment.sample_read,
        )
        sr = reads_alignment_validation.sample_read_id
        srs.add(sr)

        alignment_reference = reads_alignment.alignment_reference  # AmpliconCollectionBWAIndex
        reads_alignment.delete()
        assert not os.path.exists(reads_alignment.bam_file)
        alignment_reference.delete()
        assert not os.path.exists(alignment_reference.index_dump_dir)  # this should fail
        for f in alignment_reference.files:
            assert not os.path.exists(f)


@pytest.mark.django_db
def test_vcf_file(read_alignment):
    # parse resulting vcf and compare against test data PTBS
    for (l_id, bc_id), ra in read_alignment.items():
        # test vcf pileup
        ra_snp = pileup(ra)  # VCFReads
        assert ra_snp.reads_alignment == ra

        # test alignment of snf in vcf
        assert VCFReads.objects.filter(
            reads_alignment=ra,
        ).count() > 0
        nsp_r = VCFReads.objects.get(reads_alignment=ra)
        nsp_ra = nsp_r.reads_alignment

        nsp_r.delete()
        assert not os.path.exists(nsp_r.vcf_file)


@pytest.mark.django_db
def test_get_target_id(vcf_object):
    # parse resulting vcf and compare against test data PTBS
    for (l_id, bc_id), vcf in vcf_object.items():
        target_list = check_targets_id(vcf.vcf_file)
        assert set([6, 7]) == set(target_list)


@pytest.mark.django_db
def test_get_rel_pos(vcf_rec, snp_reads_d, requires_snps, requires_snp_targets, requires_none_genotypes):
    rel_pos = retrieve_explicit_snps_positions(vcf_rec)
    assert list(rel_pos) == [60, 65]


@pytest.mark.skipif(pytest.config.getoption("nomigrations"), reason="No migrations, no view.")
@pytest.mark.django_db
def test_parse_snps(vcf_object, snp_reads_fd, snp_reads_d, requires_snps, requires_snp_targets, requires_none_genotypes):
    # parse resulting vcf and compare against test data PTBS
    for (l_id, bc_id), vcf in vcf_object.items():
        snp_list = _parse_vcf_file(vcf.vcf_file, min_cover=1)
        for cell_snps in snp_list.keys():
            if snp_list[cell_snps]['SNP_defined'] and len(snp_list[cell_snps]['modification']) > 1:
                # print('snp_reads_fd.keys(l_id, bc_id): ', set(snp_reads_fd.keys(l_id, bc_id, int(cell_snps[0]))))
                # fs = set(snp_reads_fd.keys(l_id, bc_id, int(cell_snps[0])))
                assert snp_list[cell_snps]['modification'][0] == 'C'



@pytest.mark.skipif(pytest.config.getoption("nomigrations"), reason="No migrations, no view.")
@pytest.mark.django_db
def test_parse_vcf(vcf_object):
    for (l_id, bc_id), vcf in vcf_object.items():
        sr_snp = parse_vcf(vcf)  # SNPReads
        assert sr_snp.vcf_read == vcf

        # test alignment of snf in vcf
        assert SNPReads.objects.filter(
            vcf_read=vcf,
        ).count() > 0
        nsp_r = SNPReads.objects.get(vcf_read=vcf)
        nsp_ra = nsp_r.snps_dict
        assert nsp_r.min_cover == 1

        nsp_r.delete()
        assert not os.path.exists(nsp_r.snps_dict)

