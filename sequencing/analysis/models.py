
import os
from Bio import SeqIO
import pysam
import itertools
import re

from django.db import models
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver

from sequencing.runs.models import Demultiplexing
from targeted_enrichment.planning.models import Microsatellite, SNP
from targeted_enrichment.unwrapping.models import Unwrapper
from lib_prep.workflows.models import Library, BarcodedContent
from lib_prep.multiplexes.models import Panel


def delete_files(sender, instance, **kwargs):
    # Pass false so FileField doesn't save the model.
    for filename in getattr(instance, "files", []):
        try:
            os.unlink(filename)
        except OSError:
            raise
    for dirname in getattr(instance, "dirs", []):
        try:
            os.rmdir(dirname)
        except OSError:
            raise


class SampleReads(models.Model):
    demux = models.ForeignKey(Demultiplexing)
    barcoded_content = models.ForeignKey(BarcodedContent)
    library = models.ForeignKey(Library)
    fastq1 = models.FilePathField()
    fastq2 = models.FilePathField()

    @property
    def files(self):
        yield self.fastq1
        yield self.fastq2

post_delete.connect(delete_files, SampleReads)


class AdamMergedReads(models.Model):
    demux_reads = models.ForeignKey(SampleReads)
    # TODO: Add default path?
    # TODO: Custom FASTQ field that caches #sequences
    assembled_fastq = models.FilePathField()
    discarded_fastq = models.FilePathField()
    unassembled_forward_fastq = models.FilePathField()
    unassembled_reverse_fastq = models.FilePathField()

    @property
    def files(self):
        yield self.assembled_fastq
        yield self.discarded_fastq
        yield self.unassembled_forward_fastq
        yield self.unassembled_reverse_fastq

    def included_reads_generator(self, included_reads):
        if included_reads == 'M':
            return SeqIO.parse(self.assembled_fastq, "fastq")
        elif included_reads == 'F':
            return itertools.chain(
                SeqIO.parse(self.assembled_fastq, "fastq"),
                SeqIO.parse(self.unassembled_forward_fastq, "fastq")
            )
        else:
            raise ValueError("included_reads should be one of {}".format(AdamReadsIndex.INCLUDED_READS_OPTIONS))

post_delete.connect(delete_files, AdamMergedReads)


class AdamReadsIndex(models.Model):
    """
    A collection of reads used for generating a bowtie2 index for primers to be searched against
    """
    INCLUDED_READS_OPTIONS = (('M', 'Only merged'), ('F', 'Merged and unassembled_forward'),)
    INDEX_PREFIX = "index"
    merged_reads = models.ForeignKey(AdamMergedReads)
    included_reads = models.CharField(max_length=1, choices=INCLUDED_READS_OPTIONS)
    index_dump_dir = models.FilePathField(allow_files=False, allow_folders=True)
    padding = models.IntegerField(default=5)

    def included_reads_generator(self):
        return self.merged_reads.included_reads_generator(self.included_reads)

    @property
    def index_files_prefix(self):
        return os.path.join(self.index_dump_dir, self.INDEX_PREFIX)

    @property
    def files(self):
        for i in xrange(1, 5):
            yield "{}.{}.bt2".format(self.index_files_prefix, i)
        for i in xrange(1, 3):
            yield "{}.rev.{}.bt2".format(self.index_files_prefix, i)

    @property
    def dirs(self):
        yield self.index_dump_dir

post_delete.connect(delete_files, AdamReadsIndex)


class AdamMarginAssignment(models.Model):
    reads_index = models.ForeignKey(AdamReadsIndex)
    assignment_sam = models.FilePathField()

    def read_sam(self):
        with pysam.AlignmentFile(self.assignment_sam, "rb") as samfile:
            for r in samfile:
                if r.is_unmapped:
                    continue  # unmapped TODO: dump in appropriate bin
                read_id = samfile.getrname(r.reference_id)
                margin_name = r.query_name
                yield read_id, name_to_unwrapper_margin(margin_name)
    
    @property
    def files(self):
        yield self.assignment_sam

post_delete.connect(delete_files, AdamMarginAssignment)

LEFT, RIGHT = "left", "right"

def unwrapper_margin_to_name(unwrapper, side):
    if side not in [LEFT, RIGHT]:
        raise ValueError("Side should be one of adamiya.LEFT, adamiya.RIGHT")
    return 'unwrapper_{}_{}'.format(unwrapper.id, side)

def name_to_unwrapper_margin(name):
    m = re.match("unwrapper_([1-9][0-9]*)_({}|{})".format(LEFT, RIGHT), name)
    if not m:
        return None, None
    return Unwrapper.objects.select_subclasses().get(id=int(m.groups()[0])), m.groups()[1]


class AdamAmpliconReads(models.Model):  # This contains the actual data.
    margin_assignment = models.ForeignKey(AdamMarginAssignment)
    unwrapper = models.ForeignKey(Unwrapper)
    target_offset = models.IntegerField(null=True)
    fastq = models.FilePathField()

    class Meta:
        index_together = (
            ("margin_assignment", "unwrapper"),
        )

