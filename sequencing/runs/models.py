
from django.db import models
from lib_prep.workflows.models import Library, BarcodedContent
from primers.parts.models import IlluminaReadingAdaptor1, \
    IlluminaReadingAdaptor2


class MachineType(models.Model):
    company = models.CharField(max_length=50)
    model = models.CharField(max_length=50)
    rev_left_bc = models.BooleanField()

    def __str__(self):
        return "{}_{}".format(self.company, self.model)


class Machine(models.Model):
    machineid = models.CharField(max_length=50)
    name = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(MachineType)

    def __str__(self):
        return "{}_{}".format(self.type, self.machineid)


class NGSKit(models.Model):
    name = models.CharField(max_length=50)
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

    def __str__(self):
        return self.name


class NGSRun(models.Model):
    libraries = models.ManyToManyField(Library)
    bcl_directory = models.FilePathField(max_length=200, allow_files=False, allow_folders=True, null=True)
    name = models.CharField(max_length=100, unique=True)
    machine = models.ForeignKey(Machine)
    kit = models.ForeignKey(NGSKit, null=True)
    date = models.DateField()

    def __str__(self):
        return self.name


class DemultiplexingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()

    def __str__(self):
        return self.name


class Demultiplexing(models.Model):
    """
    A process that follows the conversion of BCL files to FASTQ by bcl2fastqe. The mulitplexed samples is identified
    by a barcode and seperates (demultiplexing).
    """
    ngs_run = models.ForeignKey(NGSRun)
    demux_scheme = models.ForeignKey(DemultiplexingScheme)

    def __str__(self):
        return "{}|{}".format(self.ngs_run, self.pk)


class MergedDemultiplexing(Demultiplexing):
    """
    A merge of few Runs to one, in order to overcome lose of data
    """
    ngs_runs = models.ManyToManyField(NGSRun)
