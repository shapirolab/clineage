
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
from targeted_enrichment.planning.models import Microsatellite, SNP
from targeted_enrichment.amplicons.models import Amplicon
from lib_prep.workflows.models import Library, BarcodedContent


def delete_files(sender, instance, **kwargs):
    # Pass false so FileField doesn't save the model.
    for filename in getattr(instance, "files", []):
        try:
            os.unlink(filename)
        except FileNotFoundError:
            pass
    for dirname in getattr(instance, "dirs", []):
        try:
            os.rmdir(dirname)
        except FileNotFoundError:
            pass


def pear_property(suffix):
    @property
    def inner(self):
        return "{}.{}.fastq".format(self.pear_files_prefix, suffix)
    return inner

class PearOutputMixin(models.Model):
    PEAR_PREFIX = "pear"
    pear_dump_dir = models.FilePathField(max_length=200, allow_files=False, allow_folders=True)

    class Meta:
        abstract = True

    @property
    def pear_files_prefix(self):
        return os.path.join(self.pear_dump_dir, self.PEAR_PREFIX)

    assembled_fastq = pear_property("assembled")
    discarded_fastq = pear_property("discarded")
    unassembled_forward_fastq = pear_property("unassembled.forward")
    unassembled_reverse_fastq = pear_property("unassembled.reverse")

    @property
    def files(self):
        yield self.assembled_fastq
        yield self.discarded_fastq
        yield self.unassembled_forward_fastq
        yield self.unassembled_reverse_fastq

    @property
    def dirs(self):
        yield self.pear_dump_dir


class BowtieIndexMixin(models.Model):
    INDEX_PREFIX = "index"
    index_dump_dir = models.FilePathField(max_length=200, allow_files=False, allow_folders=True)

    class Meta:
        abstract = True

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
    fastq1 = models.FilePathField(max_length=200)
    fastq2 = models.FilePathField(max_length=200)

    class Meta:
        index_together = (
            ("demux", "barcoded_content"),
        )
        unique_together = (
            ("demux", "barcoded_content"),
        )

    @property
    def files(self):
        yield self.fastq1
        yield self.fastq2

    def __str__(self):
        return "{} @ {}".format(self.barcoded_content.subclass, self.demux)

post_delete.connect(delete_files, SampleReads)


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
            raise ValueError("included_reads should be one of {}".format(AdamReadsIndex.INCLUDED_READS_OPTIONS))
        return itertools.filterfalse(lambda rec: re.fullmatch("N*", str(rec.seq)), it)

    def __str__(self):
        return "{}".format(self.sample_reads)

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

    def __str__(self):
        return "{}".format(self.merged_reads)

    class Meta:
        index_together = (
            ("merged_reads", "included_reads", "padding"),
        )
        unique_together = (
            ("merged_reads", "included_reads", "padding"),
        )

post_delete.connect(delete_files, AdamReadsIndex)


def _read_sam(sam_path):
    with pysam.AlignmentFile(sam_path, "rb") as samfile:
        for r in samfile:
            if r.is_unmapped:
                continue  # unmapped TODO: dump in appropriate bin
            yield samfile.getrname(r.reference_id), r.query_name


class AdamMarginAssignment(models.Model):
    reads_index = models.ForeignKey(AdamReadsIndex, unique=True)
    assignment_sam = models.FilePathField(max_length=200)

    def read_sam(self):
        for read_id, margin_name in _read_sam(self.assignment_sam):
            yield read_id, name_to_amplicon_margin(margin_name)
    
    @property
    def files(self):
        yield self.assignment_sam

    def __str__(self):
        return "{}".format(self.reads_index)

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

post_delete.connect(delete_files, AdamAmpliconReads)

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
        return "{} v. {}".format(amplicon, microsatellites_version)

post_delete.connect(delete_files, AdamMSVariations)


class Histogram(models.Model):
    sample_reads = models.ForeignKey(SampleReads)
    microsatellites_version = models.IntegerField()
    amplicon = models.ForeignKey(Amplicon)

    objects = InheritanceManager()

    # FIXME
    @property
    def subclass(self):
        return Histogram.objects.get_subclass(id=self.id)


class AdamHistogram(Histogram):
    amplicon_reads = models.ForeignKey(AdamAmpliconReads)
    assignment_sam = models.FilePathField(max_length=200)
    ms_variations = models.ForeignKey(AdamMSVariations)

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
    microsatellite = models.ForeignKey(Microsatellite, null=True)
    repeat_number = models.PositiveIntegerField()

    class Meta:
        index_together = (
            ("microsatellite", "repeat_number"),
        )
        unique_together = (
            ("microsatellite", "repeat_number"),
        )

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


class MicrosatelliteHistogramGenotypeSet(models.Model):
    microsatellite_genotype1 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')
    microsatellite_genotype2 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')
    microsatellite_genotype3 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')
    microsatellite_genotype4 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')
    microsatellite_genotype5 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')
    microsatellite_genotype6 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')
    microsatellite_genotype7 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')
    microsatellite_genotype8 = models.ForeignKey(MicrosatelliteHistogramGenotype, related_name='+')

    @staticmethod
    def genotype_field_names():
        for i in range(1,9):
            yield 'microsatellite_genotype{}'.format(i)
    
    @property
    def genotype_fields(self):
        return [
            self.microsatellite_genotype1,
            self.microsatellite_genotype2,
            self.microsatellite_genotype3,
            self.microsatellite_genotype4,
            self.microsatellite_genotype5,
            self.microsatellite_genotype6,
            self.microsatellite_genotype7,
            self.microsatellite_genotype8,
        ]

    @property
    def genotypes(self):
        for genotype in self.genotype_fields:
            if genotype.microsatellite is not None:
                yield genotype

    class Meta:
        index_together = (
            (
                "microsatellite_genotype1",
                "microsatellite_genotype2",
                "microsatellite_genotype3",
                "microsatellite_genotype4",
                "microsatellite_genotype5",
                "microsatellite_genotype6",
                "microsatellite_genotype7",
                "microsatellite_genotype8",
             ),
        )
        unique_together = (
            (
                "microsatellite_genotype1",
                "microsatellite_genotype2",
                "microsatellite_genotype3",
                "microsatellite_genotype4",
                "microsatellite_genotype5",
                "microsatellite_genotype6",
                "microsatellite_genotype7",
                "microsatellite_genotype8",
            ),
        )

class SNPHistogramGenotype(models.Model):
    snp = models.ForeignKey(SNP, null=True)
    base = models.CharField(max_length=1)

    class Meta:
        index_together = (
            ("snp", "base"),
        )
        unique_together = (
            ("snp", "base"),
        )


class SNPHistogramGenotypeSet(models.Model):
    snp_genotype1 = models.ForeignKey(SNPHistogramGenotype, related_name='+')
    snp_genotype2 = models.ForeignKey(SNPHistogramGenotype, related_name='+')
    snp_genotype3 = models.ForeignKey(SNPHistogramGenotype, related_name='+')
    snp_genotype4 = models.ForeignKey(SNPHistogramGenotype, related_name='+')
    snp_genotype5 = models.ForeignKey(SNPHistogramGenotype, related_name='+')
    snp_genotype6 = models.ForeignKey(SNPHistogramGenotype, related_name='+')
    snp_genotype7 = models.ForeignKey(SNPHistogramGenotype, related_name='+')
    snp_genotype8 = models.ForeignKey(SNPHistogramGenotype, related_name='+')

    @staticmethod
    def genotype_field_names():
        for i in range(1,9):
            yield 'snp_genotype{}'.format(i)

    @property
    def genotype_fields(self):
        return [
            self.snp_genotype1,
            self.snp_genotype2,
            self.snp_genotype3,
            self.snp_genotype4,
            self.snp_genotype5,
            self.snp_genotype6,
            self.snp_genotype7,
            self.snp_genotype8,
        ]

    @property
    def genotypes(self):
        for genotype in self.genotype_fields:
            if genotype.snp is not None:
                yield genotype

    class Meta:
        index_together = (
            (
                "snp_genotype1",
                "snp_genotype2",
                "snp_genotype3",
                "snp_genotype4",
                "snp_genotype5",
                "snp_genotype6",
                "snp_genotype7",
                "snp_genotype8",
             ),
        )
        unique_together = (
            (
                "snp_genotype1",
                "snp_genotype2",
                "snp_genotype3",
                "snp_genotype4",
                "snp_genotype5",
                "snp_genotype6",
                "snp_genotype7",
                "snp_genotype8",
             ),
        )

class HistogramEntryReads(models.Model):
    histogram = models.ForeignKey(Histogram)
    microsatellite_genotypes = models.ForeignKey(MicrosatelliteHistogramGenotypeSet)
    snp_genotypes = models.ForeignKey(SNPHistogramGenotypeSet)
    num_reads = models.PositiveIntegerField()
    fastq1 = models.FilePathField(max_length=200)
    fastq2 = models.FilePathField(max_length=200)
    fastqm = models.FilePathField(max_length=200, null=True)

    @property
    def files(self):
        yield self.fastq1
        yield self.fastq2
        yield self.fastqm

    def __str__(self):
        return "{}: {}".format(self.histogram.subclass,
            ", ".join([
                "{}".format(msg) for msg in self.microsatellite_genotypes.genotypes
            ] + [
                "{}".format(sng) for sng in self.snp_genotypes.genotypes
            ])
        )

post_delete.connect(delete_files, HistogramEntryReads)
