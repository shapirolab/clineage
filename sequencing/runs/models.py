from django.contrib.auth.models import User
from django.db import models
from lib_prep.workflows.models import Library, BarcodedContent
from primers.parts.models import IlluminaReadingAdaptor1, \
    IlluminaReadingAdaptor2

from merge import pear_with_defaults
from demux import run_bcl2fastq

class MachineType(models.Model):
    company = models.CharField(max_length=50)
    model = models.CharField(max_length=50)
    read_length = models.IntegerField(null=True)

    def __unicode__(self):
        return self.company + '_' + self.model


class Machine(models.Model):
    machineid = models.CharField(max_length=50)
    name = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(MachineType)

    def __unicode__(self):
        return self.type.__unicode__() + '_' + self.machineid


class NGSKit(models.Model):
    reading_adaptor1 = models.ForeignKey(IlluminaReadingAdaptor1)
    reading_adaptor2 = models.ForeignKey(IlluminaReadingAdaptor2)

    @property
    def adapter(self):
        """
        Return the value for the Adapter field in the SampleSheet.
        """
        return reading_adaptor2.sequence[1:].rev_comp()

    @property
    def adapterread2(self):
        """
        Return the value for the AdapterRead2 field in the SampleSheet.
        """
        return reading_adaptor1.sequence.rev_comp()


class NGSRun(models.Model):
    libraries = models.ManyToManyField(Library)
    bcl_directory = models.FilePathField(allow_files=False, allow_folders=True, null=True)
    name = models.CharField(max_length=100, unique=True)
    machine = models.ForeignKey(Machine)
    kit = models.ForeignKey(NGSKit, null=True)
    user = models.ForeignKey(User)
    date = models.DateField()


class DemultiplexingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()


class Demultiplexing(models.Model):
    ngs_run = models.ForeignKey(NGSRun)
    demux_scheme = models.ForeignKey(DemultiplexingScheme)

    def run_demux(self):
        for ss in get_sample_sheets(self):
            run_bcl2fastq(ss)


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
