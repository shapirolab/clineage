
from django.db import models
from django.contrib.contenttypes import fields
from django.contrib.auth.models import User

from genomes.models import TargetEnrichment, DNABarcode
from linapp.models import Protocol
from wet_storage.models import SampleLocation, Plate
from sampling.models import CellContent

class BarcodePair(models.Model):
    left = models.ForeignKey(DNABarcode,related_name="left")
    right = models.ForeignKey(DNABarcode,related_name="right")


class WorkFlowCell(models.Model): # cell + barcode
    content = models.ForeignKey(CellContent)
    barcodes = models.ForeignKey(BarcodePair)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')
### -------------------------------------------------------------------------------------

### -------------------------------------------------------------------------------------
class Panel(models.Model):#collection of targets  # TODO: kill and replace with PanelPlate
                                                # TODO: m2m pri_mux, well on the m2m table.
    name = models.CharField(max_length=50)
    targets = models.ManyToManyField(TargetEnrichment, related_name='panels')
    #TODO: add comment field
    def __unicode__(self):
        return self.name
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


