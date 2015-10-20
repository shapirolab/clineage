__author__ = 'ofirr'

from django.db import models
from primers.strand import StrandBaseMixin, StrandMinusMixin, StrandPlusMixin
from genomes.models import Sequence

### -------------------------------------------------------------------------------------
### Types and descriptors
### -------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------

### -------------------------------------------------------------------------------------


class KitSynthetic(models.Model,StrandBaseMixin):
    name = models.CharField(max_length=50)
    _sequence = models.ForeignKey(Sequence)

    class Meta:
        abstract = True

    @property
    def ref_sequence(self):
        return self._sequence

class ReadingAdaptorCuts(models.Model):
    overlap_start = models.IntegerField()
    overlap_end = models.IntegerField()

    @property
    def primer1tail(self):
        self.ira.sequence[self.overlap_start:]

    @property
    def overlap(self):
        self.ira.sequence[self.overlap_start:self.overlap_end]

    @property
    def primer2tail(self):
        self.ira.sequence[:self.overlap_start]

    class Meta:
        abstract = True

class IlluminaReadingAdaptor1(KitSynthetic,StrandPlusMixin):
    pass

class IlluminaReadingAdaptor2(KitSynthetic,StrandMinusMixin):
    pass

class IlluminaReadingAdaptor1Cuts(ReadingAdaptorCuts):
    ira = models.ForeignKey(IlluminaReadingAdaptor1)

class IlluminaReadingAdaptor2Cuts(ReadingAdaptorCuts):
    ira = models.ForeignKey(IlluminaReadingAdaptor2)

class IlluminaFlowCellAdaptor1(KitSynthetic,StrandPlusMixin):
    pass

class IlluminaFlowCellAdaptor2(KitSynthetic,StrandMinusMixin):
    pass

class DNABarcode1(KitSynthetic,StrandPlusMixin):
    pass

class DNABarcode2(KitSynthetic,StrandMinusMixin):
    pass
