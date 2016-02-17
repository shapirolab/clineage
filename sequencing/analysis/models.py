from django.db import models

from sequencing.runs.models import MergedReads
from targeted_enrichment.unwrapping.models import Unwrapper

from sequencing.analysis.merge import pear_with_defaults
from sequencing.analysis.demux import run_bcl2fastq


class DemultiplexingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()


class Demultiplexing(models.Model):
    ngs_run = models.ForeignKey(NGSRun)
    demux_scheme = models.ForeignKey(DemultiplexingScheme)

    def run_demux(self):
        sample_sheet_path = self.ngs_run.generate_sample_sheets()
        run_bcl2fastq(sample_sheet_path)


class DemultiplexedReads(models.Model):
    demux = models.ForeignKey(Demultiplexing)
    barcoded_content = models.ForeignKey(BarcodedContent)
    library = models.ForeignKey(Library)
    fastq1 = models.FilePathField(null=True)
    fastq2 = models.FilePathField(null=True)


class MergingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()


class MergedReads(models.Model):
    demux_read = models.ForeignKey(DemultiplexedReads)
    merge_scheme = models.ForeignKey(MergingScheme)
    # TODO: Add default path?
    # TODO: Custom FASTQ field that caches #sequences
    assembled_fastq = models.FilePathField(null=True)
    discarded_fastq = models.FilePathField(null=True)
    unassembled_forward_fastq = models.FilePathField(null=True)
    unassembled_reverse_fastq = models.FilePathField(null=True)

    def run_merge(self):
        pear_with_defaults("--forward-fastq", self.demux_read.fastq1,
                           "--reverse-fastq", self.demux_read.fastq2,
                           "--output", self.id)
        self.assembled_fastq = "{}.assembled.fastq.gz".format(self.id)
        self.discarded_fastq = "{}.discarded.fastq.gz".format(self.id)
        self.unassembled_forward_fastq = "{}.unassembled.forward.fastq.gz".format(self.id)
        self.unassembled_reverse_fastq = "{}.unassembled.reverse.fastq.gz".format(self.id)
        self.save()


class SequencingData(models.Model): # This contains the actual data.
    merged_reads = models.ForeignKey(MergedReads)
    unwrapper = models.ForeignKey(Unwrapper)
    target_offset = models.IntegerField(null=True)
    fastq = models.FilePathField(null=True)
    vcf = models.FilePathField(null=True)

    class Meta:
        index_together = (
            ("merged_reads", "unwrapper"),
        )
