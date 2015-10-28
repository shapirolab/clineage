
from django.db import models
from django.contrib.contenttypes import fields
from django.contrib.auth.models import User

from primers.models import TargetedEnrichmentReagent, PCR1PrimerPairTER, PCR1PrimerPairTERDeprecated
from linapp.models import Protocol
from wet_storage.models import SampleLocation, Plate
from sampling.models import CellContent
from primers.parts.models import DNABarcode1, DNABarcode2

class BarcodePair(models.Model):
    left = models.ForeignKey(DNABarcode1)
    right = models.ForeignKey(DNABarcode2)


class WorkFlowCell(models.Model): # cell + barcode
    content = models.ForeignKey(CellContent)
    barcodes = models.ForeignKey(BarcodePair)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')
### -------------------------------------------------------------------------------------

### -------------------------------------------------------------------------------------
### Primers Multiplex
### -------------------------------------------------------------------------------------
class TERMultiplex(models.Model): # TODO: move to primers, m2m to TER.
    name = models.CharField(max_length=20)
    primers = models.ManyToManyField(TargetedEnrichmentReagent)
    physical_locations = fields.GenericRelation(SampleLocation,
                               content_type_field='content_type',
                               object_id_field='object_id')

    class Meta:
        abstract = True

    def __unicode__(self):
        return self.name


class PCR1Multiplex(TERMultiplex):
    primers = models.ManyToManyField(PCR1PrimerPairTER)
    #TODO: physical_locations(MPXPlate)

class PCR1DeprecatedMultiplex(TERMultiplex):
    primers = models.ManyToManyField(PCR1PrimerPairTERDeprecated)
    #TODO: physical_locations(MPXPlate)

#TODO : notail?

class TargetedEnrichmentReagent(models.Model):
    te = models.ForeignKey(TargetEnrichment)
    passed_validation = models.NullBooleanField()
    validation_failure = models.ForeignKey(TargetEnrichmentFailureType, null=True)
    validation_date = models.DateField(null=True, blank=True)
    comment = models.CharField(max_length=50, blank=True, null=True)
    physical_locations = fields.GenericRelation('SampleLocation',
                                                 content_type_field='content_type',
                                                 object_id_field='object_id')

    class Meta:
        abstract = True
