
from django.db import models

from model_utils.managers import InheritanceManager

from genomes.models import DNASlice
from misc.dna import DNA
from targeted_enrichment.reagents.models import PCR1PrimerPairTER, \
    PCR1WithCompanyTagPrimerPairTER

class Amplicon(models.Model):
    slice = models.ForeignKey(DNASlice)

    objects = InheritanceManager()

    @property
    def left_margin(self):
        raise NotImplementedError()

    @property
    def right_margin(self):
        raise NotImplementedError()

class RawAmplicon(Amplicon):

    @property
    def left_margin(self):
        return DNA()

    @property
    def right_margin(self):
        return DNA()

class TargetedAmplicon(Amplicon):

    def infer_slice(self):
        raise NotImplementedError()

    class Meta:
        abstract = True

class PCR1WithCompanyTagAmplicon(TargetedAmplicon):
    ter = models.OneToOneField(PCR1WithCompanyTagPrimerPairTER)

    @property
    def left_margin(self):
        return DNA(self.ter.tag) + self.ter.te.left.ref_sequence

    @property
    def right_margin(self):
        return self.ter.te.right.ref_sequence + DNA(self.ter.tag).rev_comp()

    def infer_slice(self):
        c = self.ter.te.chromosome
        # FIXME: this is only true for 1-based,inclusive indexing!
        left = self.ter.te.left.slice.end_pos + 1
        right = self.ter.te.right.slice.start_pos - 1
        s, created = DNASlice.objects.get_or_create(
            chromosome=c,
            start_pos=left,
            end_pos=right,
        )
        self.slice = s
        #self.save()

class PCR1Amplicon(TargetedAmplicon):
    ter = models.OneToOneField(PCR1PrimerPairTER)

    @property
    def left_margin(self):
        return self.ter.te.left.ref_sequence

    @property
    def right_margin(self):
        return self.ter.te.right.ref_sequence

    def infer_slice(self):
        c = self.ter.te.chromosome
        # FIXME: this is only true for 1-based,inclusive indexing!
        left = self.ter.te.left.slice.end_pos + 1
        right = self.ter.te.right.slice.start_pos - 1
        s, created = DNASlice.objects.get_or_create(
            chromosome=c,
            start_pos=left,
            end_pos=right,
        )
        self.slice = s
        #self.save()

