
from Bio import SeqIO
import itertools
import re

from django.db import models

from sequencing.analysis.models_common import Histogram, PearOutputMixin, \
    BowtieIndexMixin, post_delete_files, SampleReads, _read_sam
from targeted_enrichment.amplicons.models import Amplicon


LEFT, RIGHT = "left", "right"

def amplicon_margin_to_name(amplicon, side):
    if side not in [LEFT, RIGHT]:
        raise ValueError("Side should be one of adamiya.LEFT, adamiya.RIGHT")
    return 'amplicon_{}_{}'.format(amplicon.id, side)
import functools

@functools.lru_cache(maxsize=100000)
def name_to_amplicon_margin(name):
    m = re.match("amplicon_([1-9][0-9]*)_({}|{})".format(LEFT, RIGHT), name)
    if not m:
        return None, None
    return Amplicon.objects.select_subclasses().get(id=int(m.groups()[0])), m.groups()[1]


class AdamMergedReads(PearOutputMixin):
    sample_reads = models.ForeignKey(SampleReads, unique=True)

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
                AdamReadsIndex.INCLUDED_READS_OPTIONS))
        return itertools.filterfalse(
            lambda rec: re.fullmatch("N*", str(rec.seq)), it
        )

    def __str__(self):
        return "{}".format(self.sample_reads)

post_delete_files(AdamMergedReads)


class AdamReadsIndex(BowtieIndexMixin):
    """
    A collection of reads used for generating a bowtie2 index for primers to be
    searched against
    """
    INCLUDED_READS_OPTIONS = (('M', 'Only merged'),
        ('F', 'Merged and unassembled_forward'),)
    merged_reads = models.ForeignKey(AdamMergedReads)
    included_reads = models.CharField(max_length=1,
        choices=INCLUDED_READS_OPTIONS)
    padding = models.IntegerField(default=5)

    def included_reads_generator(self):
        return self.merged_reads.included_reads_generator(self.included_reads)

    def __str__(self):
        return "{}".format(self.merged_reads)

    class Meta:
        index_together = (
            ("merged_reads", "included_reads", "padding"),
        )
        unique_together = (
            ("merged_reads", "included_reads", "padding"),
        )

post_delete_files(AdamReadsIndex)


class AdamMarginAssignment(models.Model):
    reads_index = models.ForeignKey(AdamReadsIndex, unique=True)
    assignment_sam = models.FilePathField(max_length=200)
    separation_finished = models.BooleanField(default=False)

    def read_sam(self):
        for read_id, margin_name in _read_sam(self.assignment_sam):
            yield read_id, name_to_amplicon_margin(margin_name)

    @property
    def files(self):
        yield self.assignment_sam

    def __str__(self):
        return "{}".format(self.reads_index)

post_delete_files(AdamMarginAssignment)


class AdamAmpliconReads(models.Model):  # This contains the actual data.
    margin_assignment = models.ForeignKey(AdamMarginAssignment)
    amplicon = models.ForeignKey(Amplicon)
    fastq1 = models.FilePathField(max_length=200)
    fastq2 = models.FilePathField(max_length=200)
    fastqm = models.FilePathField(max_length=200)

    class Meta:
        index_together = (
            ("margin_assignment", "amplicon"),
        )
        unique_together = (
            ("margin_assignment", "amplicon"),
        )

    @property
    def files(self):
        yield self.fastq1
        yield self.fastq2
        yield self.fastqm

    def __str__(self):
        return "{}[{}]".format(self.margin_assignment, self.amplicon)

post_delete_files(AdamAmpliconReads)


class AdamMSVariations(BowtieIndexMixin):
    amplicon = models.ForeignKey(Amplicon)
    padding = models.PositiveIntegerField()
    microsatellites_version = models.IntegerField()

    class Meta:
        index_together = (
            ("amplicon", "padding", "microsatellites_version"),
        )
        unique_together = (
            ("amplicon", "padding", "microsatellites_version"),
        )

    def __str__(self):
        return "{} v. {}".format(self.amplicon, self.microsatellites_version)

post_delete_files(AdamMSVariations)


class AdamHistogram(Histogram):
    amplicon_reads = models.ForeignKey(AdamAmpliconReads)
    assignment_sam = models.FilePathField(max_length=200)
    ms_variations = models.ForeignKey(AdamMSVariations)
    separation_finished = models.BooleanField(default=False)

    def read_sam(self):
        for ms_genotypes_name, read_id in _read_sam(self.assignment_sam):
            yield read_id, ms_genotypes_name

    @property
    def files(self):
        yield self.assignment_sam

    def __str__(self):
        return "{}".format(self.amplicon_reads)

    class Meta:
        index_together = (
            ("amplicon_reads", "ms_variations"),
        )
        unique_together = (
            ("amplicon_reads", "ms_variations"),
        )

post_delete_files(AdamHistogram)
