from django.contrib.auth.models import User
from django.db import models
from lib_prep.workflows.models import WorkFlowCell

__author__ = 'ofirr'

class MachineType(models.Model):
    company = models.CharField(max_length=50)
    model = models.CharField(max_length=50)

    def __unicode__(self):
        return self.company + '_' + self.model


class Machine(models.Model):
    machineid = models.CharField(max_length=50)
    name = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(MachineType)

    def __unicode__(self):
        return self.type.__unicode__() + '_' + self.machineid



class NGSRun(models.Model):
    # wfcs = models.ManyToManyField(WorkFlowCell)
    directory = models.FilePathField(null=True)
    name = models.CharField(max_length=100, unique=True)
    machine = models.ForeignKey(Machine)
    # TODO: add seqeuencing primers. kit?
    user = models.ForeignKey(User)
    date = models.DateField()


class DemultiplexingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()


class MergingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()



class DemultiplexedReads(models.Model):
    ngs_run = models.ForeignKey(NGSRun)
    directory = models.FilePathField(null=True)
    demux_scheme = models.ForeignKey(DemultiplexingScheme)


class MergedReads(models.Model):
    demux_reads = models.ForeignKey(DemultiplexedReads)
    directory = models.FilePathField(null=True)
    merge_scheme = models.ForeignKey(MergingScheme)


