__author__ = 'ofirr'

from django.db import models
from primers.strand import BaseStrandMixin, MinusStrandMixin, PlusStrandMixin
from misc.dna import DNA

### -------------------------------------------------------------------------------------
### Types and descriptors
### -------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------

### -------------------------------------------------------------------------------------

class NoStrandSynthetic(models.Model):
    name = models.CharField(max_length=50)
    _sequence = models.CharField(max_length=250)#DNAField(Sequence)

    class Meta:
        abstract = True

    @property
    def sequence(self):
        return DNA(self._sequence)

    def __str__(self):
        return "{}".format(self.name)


class KitSynthetic(models.Model, BaseStrandMixin):
    name = models.CharField(max_length=50)
    _sequence = models.CharField(max_length=250)#DNAField(Sequence)

    class Meta:
        abstract = True

    @property
    def ref_sequence(self):
        return DNA(self._sequence)

    def __str__(self):
        return "{}({})".format(self.name, self.strand)


class IlluminaReadingAdaptorForTail(models.Model, BaseStrandMixin):
    # ira = models.ForeignKey(IlluminaReadingAdaptor)
    tail_length = models.IntegerField()

    @property
    def sequence(self):
        return self.ira.sequence[-self.tail_length:]

    class Meta:
        abstract = True

    def __str__(self):
        return "{}[-{}:]".format(self.ira, self.tail_length)


class IlluminaReadingAdaptorForHead(models.Model, BaseStrandMixin):
    # ira = models.ForeignKey(IlluminaReadingAdaptor)
    head_length = models.IntegerField()

    @property
    def sequence(self):
        return self.ira.sequence[:self.head_length]

    class Meta:
        abstract = True

    def __str__(self):
        return "{}[:{}]".format(self.ira, self.head_length)


class IlluminaReadingAdaptor1(KitSynthetic, PlusStrandMixin):
    pass

class IlluminaReadingAdaptor2(KitSynthetic, MinusStrandMixin):
    pass

class IlluminaReadingAdaptor1ForTail(IlluminaReadingAdaptorForTail, PlusStrandMixin):
    ira = models.ForeignKey(IlluminaReadingAdaptor1)

class IlluminaReadingAdaptor2ForTail(IlluminaReadingAdaptorForTail, MinusStrandMixin):
    ira = models.ForeignKey(IlluminaReadingAdaptor2)

class IlluminaReadingAdaptor1ForHead(IlluminaReadingAdaptorForHead, PlusStrandMixin):
    ira = models.ForeignKey(IlluminaReadingAdaptor1)

class IlluminaReadingAdaptor2ForHead(IlluminaReadingAdaptorForHead, MinusStrandMixin):
    ira = models.ForeignKey(IlluminaReadingAdaptor2)

class IlluminaFlowCellAdaptor1(KitSynthetic, PlusStrandMixin):
    pass

class IlluminaFlowCellAdaptor2(KitSynthetic, MinusStrandMixin):
    pass

class DNABarcode1(KitSynthetic, PlusStrandMixin):
    pass

class DNABarcode2(KitSynthetic, MinusStrandMixin):
    pass

class PadlockAmplificationPlusPrimerPart1(KitSynthetic, PlusStrandMixin):
    pass

class PadlockAmplificationPlusPrimerPart2(KitSynthetic, PlusStrandMixin):
    pass

class PadlockAmplificationMinusPrimerPart1(KitSynthetic, MinusStrandMixin):
    pass

class PadlockAmplificationMinusPrimerPart2(KitSynthetic, MinusStrandMixin):
    pass

class Backbone(NoStrandSynthetic):
    pass
