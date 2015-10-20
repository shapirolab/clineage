__author__ = 'ofirr'

import re

from django.db import models
from django.contrib.auth.models import User

from genomes.models import DNASlice, Sequence, Chromosome, Target
from primers.strand import StrandBaseMixin, StrandMinusMixin, StrandPlusMixin

class UGS(models.Model,StrandBaseMixin):
    slice = models.ForeignKey(DNASlice)
    partner = models.ManyToManyField(User, null=True) # TODO: external table.

    class Meta:
        abstract = True

    @property
    def ref_sequence(self):
        return self.slice.sequence

class UGSPlus(UGS,StrandPlusMixin):
    pass

class UGSMinus(UGS,StrandMinusMixin):
    pass

class TargetEnrichment(models.Model):
    chromosome = models.ForeignKey(Chromosome)
    left = models.ForeignKey(UGSPlus)
    right = models.ForeignKey(UGSMinus)
    amplicon = models.CharField(max_length=500) # TODO: DNASlice
    targets = models.ManyToManyField(Target)
    partner = models.ManyToManyField(User, null=True) # TODO: external table.

    def update_enriched_targets(self): # return queryset of targets between the two primers and updates the m2m targets field
        assert self.left.chromosome == self.right.chromosome
        assert self.chromosome == self.left.chromosome
        self.targets = Target.objects.filter(chromosome=self.chromosome, start_pos__gte=self.left.start_pos)\
            .filter(end_pos__lte=self.right.end_pos)
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
### ----------------