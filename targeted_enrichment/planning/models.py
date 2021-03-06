__author__ = 'ofirr'

import re

from django.db import models
from django.contrib.auth.models import User

from model_utils.managers import InheritanceManager

from genomes.models import DNASlice, Chromosome, RestrictionSiteDNASlice
from misc.dna import DNA
from primers.strand import BaseStrandMixin, MinusStrandMixin, PlusStrandMixin


class UGS(models.Model,BaseStrandMixin):
    slice = models.ForeignKey(DNASlice, unique=True)

    class Meta:
        abstract = True

    @property
    def ref_sequence(self):
        return self.slice.sequence

    def __str__(self):
        return "{}({})".format(self.slice, self.strand)

    def __len__(self):
        return len(self.slice)


class UGSPlus(UGS, PlusStrandMixin):
    pass


class UGSMinus(UGS, MinusStrandMixin):
    pass


class Target(models.Model):
    name = models.CharField(max_length=50)
    slice = models.ForeignKey(DNASlice)
    partner = models.ManyToManyField(User) # TODO: external table.
    
    objects = InheritanceManager

    def __str__(self):
        return "{} @ {}".format(self.name, self.slice)


class TargetEnrichment(models.Model):
    chromosome = models.ForeignKey(Chromosome, null=True)
    left = models.ForeignKey(UGSPlus)
    right = models.ForeignKey(UGSMinus)
    # slice = models.ForeignKey(DNASlice)
    partner = models.ManyToManyField(User) # TODO: external table.
    planning_version = models.IntegerField()

    def __str__(self):
        return "{}, {}".format(self.left, self.right)


class RestrictionEnzyme(models.Model):  # repopulate from scratch, no migration
    name = models.CharField(max_length=50)
    _sequence = models.CharField(max_length=50) # TODO: DNAField
    cut_delta = models.IntegerField()  # position of cutting site relative to start_pos
    sticky_bases = models.IntegerField()
    sequence_len = models.PositiveIntegerField()

    @property
    def sequence(self):
        return DNA(self._sequence)

    def save(self, *args, **kwargs):
        self.sequence_len = len(self.sequence)
        return super(RestrictionEnzyme, self).save(*args, **kwargs)

    def __str__(self):
        return self.name


class RestrictionSite(models.Model):
    slice = models.ForeignKey(RestrictionSiteDNASlice, db_column='slice_id')
    enzyme = models.ForeignKey(RestrictionEnzyme, related_name="sites")

    @property
    def sequence(self):
        return self.enzyme.sequence

    def __str__(self):
        return "{} @ {}".format(self.enzyme.name, self.slice)


class Microsatellite(Target):
    repeat_unit_len = models.PositiveIntegerField() #length of repeat Nmer
    repeat_unit_type = models.CharField(max_length=50) #string of repeat Nmer
    repeat_number = models.DecimalField(max_digits=5, decimal_places=1, null=True)
    repeat_unit_ref_seq = models.CharField(max_length=50) #string of acutal unit in genome.
    planning_version = models.IntegerField()

    def __str__(self):
        return "{}x{} @ {}".format(self.repeat_number, self.repeat_unit_type,
            self.slice)


class SNP(Target):
    mutation = models.CharField(max_length=10, null=True) #X>Y
    modified = models.CharField(max_length=10, null=True) #Y

    def __str__(self):
        return "{}:{} @ {}".format(self.name, self.mutation, self.slice)


# TODO: add indel
