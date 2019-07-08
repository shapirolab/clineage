
from Bio import SeqIO
import itertools
import re

from django.db import models

from lib_prep.multiplexes.models import OM6Panel
from sequencing.analysis.models_common import Histogram, PearOutputMixin, \
    BowtieIndexMixin, post_delete_files, SampleReads, _read_bam
from targeted_enrichment.amplicons.models import Amplicon, AmpliconCollection


class FullMSVMergedReads(PearOutputMixin):
    sample_reads = models.ForeignKey(SampleReads, unique=True)
    _num_reads_M = models.PositiveIntegerField(null=True, default=None)
    _num_reads_F = models.PositiveIntegerField(null=True, default=None)
    INCLUDED_READS_OPTIONS = (('M', 'Only merged'),
                              ('F', 'Merged and unassembled_forward'),)


    def num_reads(self, included_reads):
        attr_name = '_num_reads_{}'.format(included_reads)
        if included_reads not in ['M', 'F']:
            raise ValueError("included_reads should be one of {}".format(
                FullMSVMergedReads.INCLUDED_READS_OPTIONS))
        if getattr(self, attr_name) is None:
            num_reads = sum(1 for x in self.included_reads_generator(included_reads))
            setattr(self,attr_name, num_reads)
            self.save()
        return getattr(self, attr_name)


    def included_reads_generator(self, included_reads):
        if included_reads == 'M':
            it = SeqIO.parse(self.assembled_fastq, "fastq")
        elif included_reads == 'F':
            it = itertools.chain(
                SeqIO.parse(self.assembled_fastq, "fastq"),
                SeqIO.parse(self.unassembled_forward_fastq, "fastq")
            )
        else:
            raise ValueError("included_reads should be one of {}".format(
                FullMSVMergedReads.INCLUDED_READS_OPTIONS))
        return itertools.filterfalse(
            lambda rec: re.fullmatch("N*", str(rec.seq)), it
        )

    def __str__(self):
        return "{}".format(self.sample_reads)

post_delete_files(FullMSVMergedReads)


class FullMSVMergedReadsPart(models.Model):
    merged_reads = models.ForeignKey(FullMSVMergedReads)
    fastq_part = models.FilePathField(max_length=200)
    start_row = models.IntegerField()
    rows = models.IntegerField()

    @property
    def files(self):
        yield self.fastq_part

    def __str__(self):
        return "{} part {}+{}".format(self.merged_reads, self.start_row, self.rows)

    class Meta:
        unique_together = (
            ("merged_reads", "start_row", "rows"),
        )
post_delete_files(FullMSVMergedReadsPart)


class FullMSVariations(BowtieIndexMixin):
    amplicon_collection = models.ForeignKey(AmpliconCollection)
    padding = models.PositiveIntegerField()
    microsatellites_version = models.IntegerField()

    class Meta:
        index_together = (
            ("amplicon_collection", "padding", "microsatellites_version"),
        )
        unique_together = (
            ("amplicon_collection", "padding", "microsatellites_version"),
        )

    def __str__(self):
        return "{} v. {}".format(self.amplicon_collection, self.microsatellites_version)

post_delete_files(FullMSVariations)


class FullMSVAssignmentPart(models.Model):
    merged_reads_part = models.ForeignKey(FullMSVMergedReadsPart)
    assignment_bam = models.FilePathField(max_length=200)
    ms_variations = models.ForeignKey(FullMSVariations)

    class Meta:
        unique_together = (
            ("merged_reads_part", "ms_variations"),
        )

    def read_bam(self):
        for ms_genotypes_name, read_id in _read_bam(self.assignment_bam):
            yield read_id, ms_genotypes_name

    @property
    def files(self):
        yield self.assignment_bam

    def __str__(self):
        return "{}@{}".format(self.merged_reads_part, self.ms_variations.amplicon_collection)

post_delete_files(FullMSVAssignmentPart)


class FullMSVAssignment(models.Model):
    merged_reads = models.ForeignKey(FullMSVMergedReads)
    sorted_assignment_bam = models.FilePathField(max_length=200)
    ms_variations = models.ForeignKey(FullMSVariations)
    separation_finished = models.BooleanField(default=False)

    class Meta:
        unique_together = (
            ("merged_reads", "ms_variations"),
        )

    def read_bam(self):
        for ms_genotypes_name, read_id in _read_bam(self.sorted_assignment_bam):
            yield read_id, ms_genotypes_name

    @property
    def files(self):
        yield self.sorted_assignment_bam

    def __str__(self):
        return "{}@{}".format(self.merged_reads, self.ms_variations.amplicon_collection)

post_delete_files(FullMSVAssignment)


class FullMSVHistogram(Histogram):
    assignment = models.ForeignKey(FullMSVAssignment)
    amplicon_copy = models.ForeignKey(Amplicon)  # to allow unique and index together, amplicon field and assignment field must reside on the same table

    def __str__(self):
        return "{}".format(self.assignment)

    class Meta:
        index_together = (
            ("amplicon_copy", "assignment"),
        )
        unique_together = (
            ("amplicon_copy", "assignment"),
        )

post_delete_files(FullMSVHistogram)
