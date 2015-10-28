
from django.db import models
from django.contrib.contenttypes import fields
from django.contrib.auth.models import User

from genomes.models import TargetEnrichment
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
class PrimersMultiplex(models.Model): # TODO: move to primers, m2m to TER.
    name = models.CharField(max_length=20)
    primers = models.ManyToManyField(TargetEnrichment)
    physical_locations = fields.GenericRelation(SampleLocation,
                               content_type_field='content_type',
                               object_id_field='object_id')

    def __unicode__(self):
        return self.name


