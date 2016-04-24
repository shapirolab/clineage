
import os
from Bio import SeqIO
import pysam
import itertools
import re

from model_utils.managers import InheritanceManager

from django.db import models
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver

from sequencing.runs.models import Demultiplexing
from targeted_enrichment.planning.models import Microsatellite, SNP, \
    PhasedMicrosatellites
from targeted_enrichment.amplicons.models import Amplicon
from lib_prep.workflows.models import Library, BarcodedContent


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


class BowtieIndexMixin(models.Model):
    INDEX_PREFIX = "index"
    index_dump_dir = models.FilePathField(allow_files=False, allow_folders=True)

    class Meta:
        abstract=True

    @property
    def index_files_prefix(self):
        return os.path.join(self.index_dump_dir, self.INDEX_PREFIX)

    @property
    def files(self):
        for i in range(1, 5):
            yield "{}.{}.bt2".format(self.index_files_prefix, i)
        for i in range(1, 3):
            yield "{}.rev.{}.bt2".format(self.index_files_prefix, i)

    @property
    def dirs(self):
        yield self.index_dump_dir


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
    sample_reads = models.ForeignKey(SampleReads)
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


class AdamReadsIndex(BowtieIndexMixin):
    """
    A collection of reads used for generating a bowtie2 index for primers to be searched against
    """
    INCLUDED_READS_OPTIONS = (('M', 'Only merged'), ('F', 'Merged and unassembled_forward'),)
    merged_reads = models.ForeignKey(AdamMergedReads)
    included_reads = models.CharField(max_length=1, choices=INCLUDED_READS_OPTIONS)
    padding = models.IntegerField(default=5)

    def included_reads_generator(self):
        return self.merged_reads.included_reads_generator(self.included_reads)

post_delete.connect(delete_files, AdamReadsIndex)


def _read_sam(sam_path):
    with pysam.AlignmentFile(sam_path, "rb") as samfile:
        for r in samfile:
            if r.is_unmapped:
                continue  # unmapped TODO: dump in appropriate bin
            yield samfile.getrname(r.reference_id), r.query_name


class AdamMarginAssignment(models.Model):
    reads_index = models.ForeignKey(AdamReadsIndex)
    assignment_sam = models.FilePathField()

    def read_sam(self):
        for read_id, margin_name in _read_sam(self.assignment_sam):
            yield read_id, name_to_amplicon_margin(margin_name)
    
    @property
    def files(self):
        yield self.assignment_sam

post_delete.connect(delete_files, AdamMarginAssignment)


LEFT, RIGHT = "left", "right"

def amplicon_margin_to_name(amplicon, side):
    if side not in [LEFT, RIGHT]:
        raise ValueError("Side should be one of adamiya.LEFT, adamiya.RIGHT")
    return 'amplicon_{}_{}'.format(amplicon.id, side)

def name_to_amplicon_margin(name):
    m = re.match("amplicon_([1-9][0-9]*)_({}|{})".format(LEFT, RIGHT), name)
    if not m:
        return None, None
    return Amplicon.objects.select_subclasses().get(id=int(m.groups()[0])), m.groups()[1]


class AdamAmpliconReads(models.Model):  # This contains the actual data.
    margin_assignment = models.ForeignKey(AdamMarginAssignment)
    amplicon = models.ForeignKey(Amplicon)
    target_offset = models.IntegerField(null=True)
    fastq1 = models.FilePathField()
    fastq2 = models.FilePathField()
    fastqm = models.FilePathField()

    class Meta:
        index_together = (
            ("margin_assignment", "amplicon"),
        )

    @property
    def files(self):
        yield self.fastq1
        yield self.fastq2
        yield self.fastqm

post_delete.connect(delete_files, AdamAmpliconReads)

class AdamMSVariations(BowtieIndexMixin):
    amplicon = models.ForeignKey(Amplicon)
    padding = models.PositiveIntegerField()
    microsatellites_version = models.ForeignKey(PhasedMicrosatellites)

    class Meta:
        index_together=[
            ("amplicon","padding")
        ]

post_delete.connect(delete_files, AdamMSVariations)


class Histogram(models.Model):
    sample_reads = models.ForeignKey(SampleReads)

    objects = InheritanceManager()


class AdamHistogram(Histogram):
    amplicon_reads = models.ForeignKey(AdamAmpliconReads)
    assignment_sam = models.FilePathField()
    ms_variations = models.ForeignKey(AdamMSVariations)

    def read_sam(self):
        for ms_genotypes_name, read_id in _read_sam(self.assignment_sam):
            ms_genotypes, prefix = name_to_ms_genotypes(ms_genotypes_name)
            assert int(prefix) == self.amplicon_reads.amplicon_id
            yield read_id, ms_genotypes
    
    @property
    def files(self):
        yield self.assignment_sam

post_delete.connect(delete_files, AdamHistogram)


def ms_genotypes_to_name(ms_genotypes, prefix):
    names = ["{}".format(mhg) for mhg in ms_genotypes]
    return ":".join([prefix] + names)


def name_to_ms_genotypes(ms_genotypes_name):
    genotypes_plus = ms_genotypes_name.split(":")
    prefix = genotypes_plus[0]
    ms_genotypes = tuple([MicrosatelliteHistogramGenotype.get_for_string(s) \
        for s in genotypes_plus[1:]])
    return ms_genotypes, prefix


class MicrosatelliteHistogramGenotype(models.Model):
    microsatellite = models.ForeignKey(Microsatellite)
    repeat_number = models.PositiveIntegerField()

    def __str__(self):
        return "{}={}".format(self.microsatellite.id, self.repeat_number)

    @classmethod
    def get_for_string(cls, s):
        # py3: m = re.fullmatch("([1-9][0-9]*)=([1-9][0-9]*)", s)
        m = re.match("([1-9][0-9]*)=([1-9][0-9]*)$", s)
        if not m:
            raise ValueError("Bad ms genotype string: {}".format(s))
        msid, rn = m.groups()
        # NOTE: we save a query by not getting the actual MS.
        # This is OK as long as we use a db with FK enforcement.
        obj, c = cls.objects.get_or_create(
            microsatellite_id=int(msid),
            repeat_number=int(rn),
        )
        return obj

    @classmethod
    def get_for_genotype(cls, ms, rn):
        mhg, c = cls.objects.get_or_create(
            microsatellite=ms,
            repeat_number=rn,
        )
        return mhg

    @property
    def sequence(self):
        return self.microsatellite.repeat_unit_ref_seq * self.repeat_number


class SNPHistogramGenotype(models.Model):
    snp = models.ForeignKey(SNP)
    base = models.CharField(max_length=1)


class HistogramEntryReads(models.Model):
    histogram = models.ForeignKey(Histogram)
    amplicon = models.ForeignKey(Amplicon)
    microsatellites_version = models.ForeignKey(PhasedMicrosatellites)
    microsatellite_genotypes = models.ManyToManyField(MicrosatelliteHistogramGenotype)
    snp_genotypes = models.ManyToManyField(SNPHistogramGenotype)
    num_reads = models.PositiveIntegerField()
    fastq1 = models.FilePathField()
    fastq2 = models.FilePathField()
    fastqm = models.FilePathField(null=True)

    @property
    def files(self):
        yield self.fastq1
        yield self.fastq2
        yield self.fastqm

post_delete.connect(delete_files, HistogramEntryReads)
