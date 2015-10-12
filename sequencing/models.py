import os

from django.db import models
from django.conf import settings
from django.dispatch import receiver
from django.db.models.signals import post_save
from django.contrib.auth.models import User

from sampling.models import CellContent
from lib_prep.models import WorkFlowCell
from genomes.models import Target

### -------------------------------------------------------------------------------------

class NGSRun(models.Model):
    wfcs = models.ManyToManyField(WorkFlowCell)
    directory = models.FilePathField(null=True)
    name = models.CharField(max_length=100, unique=True)
    machine = models.ForeignKey(Machine)
    # TODO: add seqeuencing primers. kit?
    user = models.ForeignKey(User)
    date = models.DateField()

class DemultiplexedReads(models.Model):
    ngs_run = models.ForeignKey(NGSRun)
    directory = models.FilePathField(null=True)
    demux_scheme = models.ForeignKey(DemultiplexingScheme)

class MergedReads(models.Model):
    demux_reads = models.ForeignKey(DemultiplexedReads)
    directory = models.FilePathField(null=True)
    merge_scheme = models.ForeignKey(MergingScheme)
### -------------------------------------------------------------------------------------
class DemultiplexingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()

class MergingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()

class MachineType(models.Model):
    company = models.CharField(max_length=50)
    model = models.CharField(max_length=50)

    def __unicode__(self):
        return self.company + '_' + self.model

class Machine(models.Model):
    machineid = models.CharField(max_length=50)
    name = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(MachineType)
    ip = models.IPAddressField(null=True, blank=True)

    def __unicode__(self):
        return self.type.__unicode__() + '_' + self.machineid

### -------------------------------------------------------------------------------------
class SequencingData(models.Model): # This contains the actual data.
    cell_content = models.ForeignKey(CellContent)
    merged_reads = models.ForeignKey(MergedReads)
    target = models.ForeignKey(Target)
    target_offset = models.IntegerField(null=True)
    fastq = models.FilePathField(null=True)
    vcf = models.FilePathField(null=True)

    class Meta:
        index_together = (
            ("cell_content", "merged_reads", "target")
        )
