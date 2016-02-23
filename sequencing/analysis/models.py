from django.db import models

from sequencing.runs.models import NGSRun
from targeted_enrichment.unwrapping.models import Unwrapper
from lib_prep.workflows.models import Library, BarcodedContent

from sequencing.analysis.merge import pear_with_defaults
from sequencing.analysis.demux import run_bcl2fastq

from Bio import SeqIO
from itertools import chain


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

    @property
    def pear_prefix(self):
        return self.id

    def run_merge(self):
        pear_with_defaults("-f", self.demux_read.fastq1,
                           "-r", self.demux_read.fastq2,
                           "-o", self.pear_prefix)
        self.assembled_fastq = "{}.assembled.fastq".format(self.pear_prefix)
        self.discarded_fastq = "{}.discarded.fastq".format(self.pear_prefix)
        self.unassembled_forward_fastq = "{}.unassembled.forward.fastq".format(self.pear_prefix)
        self.unassembled_reverse_fastq = "{}.unassembled.reverse.fastq".format(self.pear_prefix)
        self.save()


class ReadsIndex(models.Model):
    """
    A collection of reads used for generating a bowtie2 index for primers to be searched against
    """
    INCLUDED_READS_OPTIONS = (('M', 'Only merged'), ('F', 'Merged and unassembled_forward'),)
    merged_reads = models.ForeignKey(MergedReads)
    included_reads = models.CharField(max_length=1, choices=INCLUDED_READS_OPTIONS)
    padded_reads_fasta = models.FilePathField(null=True)
    padding = models.IntegerField(default=5)

    @property
    def padded_reads_fasta_name(self):
        return "{}.padded.fa".format(self.merged_reads.pear_prefix)

    def pad_records(self, records):
        for record in records:
            yield "N"*self.padding + record + "N"*self.padding

    def included_reads_generator(self):
        if self.included_reads == 'M':
            return SeqIO.parse(self.merged_reads.assembled_fastq, "fastq")
        if self.included_reads == 'F':
            return chain(
                SeqIO.parse(self.merged_reads.assembled_fastq, "fastq"),
                SeqIO.parse(self.merged_reads.unassembled_forward_fastq, "fastq")
            )
        #IntegrityError
        raise Exception()

    def create_final_merged_fastq(self):
        padded_reads = self.pad_records(self.included_reads_generator())
        SeqIO.write(padded_reads, self.padded_reads_fasta_name, "fasta")
        self.padded_reads_fasta = self.padded_reads_fasta_name
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
