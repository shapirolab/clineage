__author__ = 'ofirr'

import re

from django.db import models

from django.contrib.auth.models import User

from genomes.models import DNASlice, Chromosome
from primers.strand import BaseStrandMixin, MinusStrandMixin, PlusStrandMixin


class UGS(models.Model,BaseStrandMixin):
    slice = models.ForeignKey(DNASlice)
    partner = models.ManyToManyField(User, null=True) # TODO: external table.

    class Meta:
        abstract = True

    @property
    def ref_sequence(self):
        return self.slice.sequence

class UGSPlus(UGS,PlusStrandMixin):
    pass

class UGSMinus(UGS,MinusStrandMixin):
    pass

class Target(models.Model):
    name = models.CharField(max_length=50)
    slice = models.ForeignKey(DNASlice)
    partner = models.ManyToManyField(User, null=True) # TODO: external table.

class TargetEnrichment(models.Model):
    chromosome = models.ForeignKey(Chromosome)
    left = models.ForeignKey(UGSPlus)
    right = models.ForeignKey(UGSMinus)
    amplicon = models.CharField(max_length=500) # TODO: unwrapper/DNASlice
    # slice = models.ForeignKey(DNASlice)
    targets = models.ManyToManyField(Target)
    partner = models.ManyToManyField(User, null=True) # TODO: external table.

    def update_enriched_targets(self): # return queryset of targets between the two primers and updates the m2m targets field
        assert self.left.slice.chromosome == self.right.slice.chromosome
        assert self.chromosome == self.left.slice.chromosome
        # TODO: change to slice query.
        self.targets = Target.objects.filter(slice__chromosome=self.chromosome, slice__start_pos__gte=self.left.slice.start_pos)\
            .filter(slice__end_pos__lte=self.right.slice.end_pos)
        self.save()
        return self.targets.all()

    def amplicon_indices(self): # TODO: kill
        return (self.left.start_pos, self.right.end_pos)

    def get_internal_restriction(self, restriction):
        return [self.amplicon_indices()[0] + m.start() for m in re.finditer(restriction, self.chromosome.getdna(*self.amplicon_indices()))]

    def get_surrounding_restriction(self, restriction, max_seek=5000):
        for x in range(0, max_seek, 10):
            lamplicon = self.chromosome.getdna(self.amplicon_indices()[0]-x, self.amplicon_indices()[0])
            lttaas = [self.amplicon_indices()[0] - m.start() for m in re.finditer(restriction, lamplicon)]
            if lttaas:
                break

        for x in range(0, max_seek, 10):
            ramplicon = self.chromosome.getdna(self.amplicon_indices()[1], self.amplicon_indices()[1]+x)
            rttaas = [self.amplicon_indices()[1] + m.start() for m in re.finditer(restriction, ramplicon)]
            if rttaas:
                break

        if lttaas and rttaas:
            return max(lttaas), min(rttaas)
        return None


    def __unicode__(self):
        return 'TE: left=%s, right=%s' % (self.left.name, self.right.name)

class RestrictionEnzyme(models.Model):  # repopulate from scratch, no migration
    name = models.CharField(max_length=50)
    sequence = models.CharField(max_length=50) # TODO: DNAField
    cut_delta = models.IntegerField()  # position of cutting site relative to start_pos
    sticky_bases = models.IntegerField()
    sequence_len = models.PositiveIntegerField()

    def save(self, *args, **kwargs):
        self.sequence_len = len(self.sequence)
        return super(RestrictionEnzyme, self).save(*args, **kwargs)

    def __unicode__(self):
        return self.name

class RestrictionSite(models.Model):
    slice = models.ForeignKey(DNASlice)
    enzyme = models.ForeignKey(RestrictionEnzyme, related_name="sites")

    @property
    def sequence(self):
        return self.enzyme.sequence

class Microsatellite(Target):
    repeat_unit_len = models.PositiveIntegerField() #length of repeat Nmer
    repeat_unit_type = models.CharField(max_length=50) #string of repeat Nmer
    repeat_number = models.DecimalField(max_digits=5, decimal_places=1, null=True)


class SNP(Target):
    mutation = models.CharField(max_length=10) #X>Y
    modified = models.CharField(max_length=10) #Y

# TODO: add indel
