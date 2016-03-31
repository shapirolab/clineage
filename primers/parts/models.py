__author__ = 'ofirr'

from django.db import models
from primers.strand import BaseStrandMixin, MinusStrandMixin, PlusStrandMixin
from misc.dna import DNA

### -------------------------------------------------------------------------------------
### Types and descriptors
### -------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------

### -------------------------------------------------------------------------------------


class KitSynthetic(models.Model,BaseStrandMixin):
    name = models.CharField(max_length=50)
    _sequence = models.CharField(max_length=250)#DNAField(Sequence)

    class Meta:
        abstract = True

    @property
    def ref_sequence(self):
        return DNA(self._sequence)

    def __unicode__(self):
        return "{}({})".format(self.name, self.strand)

class ReadingAdaptorCuts(models.Model):
    overlap_start = models.IntegerField()
    overlap_end = models.IntegerField()

    @property
    def primer1tail(self):
        return self.ira.sequence[self.overlap_start:]

    @property
    def overlap(self):
        return self.ira.sequence[self.overlap_start:self.overlap_end]

    @property
    def primer2tail(self):
        return self.ira.sequence[:self.overlap_start]

    class Meta:
        abstract = True

    def __unicode__(self):
        return "{}[:{}:{}:]".format(self.ira, self.overlap_start,
            self.overlap_end)

class IlluminaReadingAdaptor1(KitSynthetic,PlusStrandMixin):
    pass

class IlluminaReadingAdaptor2(KitSynthetic,MinusStrandMixin):
    pass

class IlluminaReadingAdaptor1Cuts(ReadingAdaptorCuts):
    ira = models.ForeignKey(IlluminaReadingAdaptor1)

class IlluminaReadingAdaptor2Cuts(ReadingAdaptorCuts):
    ira = models.ForeignKey(IlluminaReadingAdaptor2)

class IlluminaFlowCellAdaptor1(KitSynthetic,PlusStrandMixin):
    pass

class IlluminaFlowCellAdaptor2(KitSynthetic,MinusStrandMixin):
    pass

class DNABarcode1(KitSynthetic,PlusStrandMixin):
    pass

class DNABarcode2(KitSynthetic,MinusStrandMixin):
    pass

class PadlockAmplificationPlusPrimer(KitSynthetic,PlusStrandMixin):
    pass

class PadlockAmplificationMinusPrimer(KitSynthetic,MinusStrandMixin):
    pass
