
from django.contrib.auth.models import User
from django.db import models
from django.db.models.query import QuerySet
from lib_prep.workflows.models import Library, BarcodedContent
from primers.parts.models import IlluminaReadingAdaptor1, \
    IlluminaReadingAdaptor2

from sequencing.runs.demux import run_bcl2fastq, generate_sample_sheets


class MachineType(models.Model):
    company = models.CharField(max_length=50)
    model = models.CharField(max_length=50)

    def __str__(self):
        return self.company + '_' + self.model


class Machine(models.Model):
    machineid = models.CharField(max_length=50)
    name = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(MachineType)

    def __str__(self):
        return self.type.__str__() + '_' + self.machineid


class NGSKit(models.Model):
    reading_adaptor1 = models.ForeignKey(IlluminaReadingAdaptor1)
    reading_adaptor2 = models.ForeignKey(IlluminaReadingAdaptor2)
    read_length = models.IntegerField(null=True)

    @property
    def fwd_read_adaptor(self):
        """
        Return the value for the Adapter field in the SampleSheet.
        """
        return self.reading_adaptor2.sequence[1:].rev_comp()

    @property
    def rev_read_adaptor(self):
        """
        Return the value for the AdapterRead2 field in the SampleSheet.
        """
        return self.reading_adaptor1.sequence.rev_comp()


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
