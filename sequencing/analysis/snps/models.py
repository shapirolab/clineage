# __author__ = 'veronika'

import pickle
import vcf
from Bio.SeqRecord import SeqRecord

from django.db import models
from sequencing.analysis.snps.parse_snps import retrieve_explicit_snps_positions, snp_is_ms
from targeted_enrichment.amplicons.models import Amplicon, AmpliconCollection, PlainTargetedAmplicon
from sequencing.analysis.models_common import SampleReads, _read_bam, post_delete_files, BWAIndexMixin


def is_snp(r, min_cover=5):
    if r.INFO['DP'] < min_cover or len(r.REF) != 1:
        return False
    # if len(r.ALT) == 1:
    #         return False
    return True


def _parse_vcf_file(vcf_file, min_cover=5):
    cell_snps = {}
    chrom = None
    is_ms = False
    snp_rel_pos = []

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for r in vcf_reader:
        if is_snp(r, min_cover):
            if chrom != r.CHROM:
                chrom = r.CHROM
                snp_rel_pos = list(retrieve_explicit_snps_positions(r))
                is_ms = snp_is_ms(r)
            if is_ms:
                continue
            smp_mod = r.ALT
            if len(r.ALT) == 1:
                smp_mod = SeqRecord('X')
            cell_snps.setdefault((r.CHROM, str(r.POS)), {'base': r.REF,
                                                         'modification': smp_mod,
                                                         'DP': r.INFO['DP'],
                                                         'stats': r.INFO['QS'],
                                                         'SNP_defined': False})

            for snp_rel in snp_rel_pos:
                if r.POS == snp_rel:
                    cell_snps[(r.CHROM, str(r.POS))]['SNP_defined'] = True
                    # Note that vcf is 1-based, so we should have abs_pos-amplicon_start+1 as rel_pos.
                    # However, because bedtools is 0-based, our amplicons are without the first base,
                    # so vcf indices become 0-based in the original amplicon, and rel_pos = abs_pos-amplicon_start.
    return cell_snps


def _extract_snp_file(snp_file):
    with open(snp_file, 'rb') as handle:
        snp_dict = pickle.load(handle)
    return snp_dict


class AmpliconCollectionBWAIndex(BWAIndexMixin):
    amplicon_collection = models.ForeignKey(AmpliconCollection, unique=True)

    def __str__(self):
        return "{}".format(self.amplicon_collection)

post_delete_files(AmpliconCollectionBWAIndex)


class ReadsAlignment(models.Model):
    sample_read = models.ForeignKey(SampleReads)
    alignment_reference = models.ForeignKey(AmpliconCollectionBWAIndex)
    bam_file = models.FilePathField(max_length=200, allow_files=True, allow_folders=False)

    def parse_bam_file(self):
        for ms_genotypes_name, read_id in _read_bam(self.bam_file):
            yield read_id, ms_genotypes_name

    @property
    def files(self):
        yield self.bam_file

    def __str__(self):
        return "{}@{}".format(self.sample_read, self.alignment_reference)

    class Meta:
        unique_together = (
            ("sample_read", "alignment_reference"),
        )
post_delete_files(ReadsAlignment)


class VCFReads(models.Model):
    reads_alignment = models.ForeignKey(ReadsAlignment, unique=True)
    vcf_file = models.FilePathField(max_length=200, allow_files=True, allow_folders=False)

    def parse_vcf_file(self):
        rows = _parse_vcf_file(self.vcf_file)
        return rows

    @property
    def files(self):
        yield self.vcf_file

    def __str__(self):
        return "{}@{}".format(self.reads_alignment, self.vcf_file)

post_delete_files(VCFReads)


class SNPReads(models.Model):
    vcf_read = models.ForeignKey(VCFReads)
    min_cover = models.IntegerField()
    snps_dict = models.FilePathField(max_length=200, allow_files=True, allow_folders=False)

    def extract_snp_file(self):
        snp_dict = _extract_snp_file(self.snps_dict)
        return snp_dict

    @property
    def files(self):
        yield self.snps_dict

    def __str__(self):
        return "{}@{}".format(self.vcf_read, self.snps_dict)

    class Meta:
        unique_together = (
            ("vcf_read", "min_cover"),
        )
post_delete_files(SNPReads)
