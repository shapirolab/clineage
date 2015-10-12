
from django.db import models
from django.contrib.contenttypes import fields
from django.contrib.auth.models import User

from genomes.models import TargetEnrichment, DNABarcode
from linapp.models import Protocol
from wet_storage.models import SampleLocation, Plate
from sampling.models import CellContent

class WorkFlowCell(models.Model): # cell + barcode
    content = models.ForeignKey(CellContent)
    barcodes = models.ForeignKey(BarcodePair)
    physical_locations = fields.GenericRelation('SampleLocation',
                                             content_type_field='content_type',
                                             object_id_field='object_id')

class TamirLibCells(models.Model):
    wfc = models.ForeignKey(TamirWFC)
    lib = models.ForeignKey(TamirLib)
    volume = models.FloatField()

class TamirWFC(WorkFlowCell):
    workflow48 = models.ForeignKey(TamirWF48)
    aar_well = models.CharField(max_length=3)
    post_pcr1_conc = models.FloatField()
    post_pcr1_location = fields.GenericRelation('SampleLocation',
                                             content_type_field='content_type',
                                             object_id_field='object_id')
    pre_pcr2_conc = models.FloatField()
    pre_pcr2_location = fields.GenericRelation('SampleLocation',
                                             content_type_field='content_type',
                                             object_id_field='object_id')
    barcode_well = models.CharField(max_length=3)
    dilution_well = models.CharField(max_length=3)
    final_well = models.CharField(max_length=3)

    class Meta:
        unique_together = (
            ("workflow48", "aar_well"),
            ("workflow48", "barcode_well"),
            ("workflow48", "dilution_well"),
            ("workflow48", "final_well"),
            # can't do uniqueness with the larger workflow
        )
### -------------------------------------------------------------------------------------
class TamirWF48(models.Model):
    workflow384 = models.ForeignKey(TamirWF384)
    #TODO: enforce plate types
    aar = models.ForeignKey(Plate)
    # TODO: panel?
    barcode_plate = models.ForeignKey(Plate)
    dilution_plate = models.ForeignKey(Plate)

class TamirWF384(models.Model):
    #TODO: enforce plate type
    final_plate = models.ForeignKey(Plate)
    bluepippin_length = models.IntegerField()

### -------------------------------------------------------------------------------------
class BarcodePair(models.Model):
    left = models.ForeignKey(DNABarcode)
    right = models.ForeignKey(DNABarcode)
### -------------------------------------------------------------------------------------
class BarcodePlate(models.Model):
    #TODO: enforce plate type
    plate = models.ForeignKey(Plate)
    barcodes = models.ManyToManyField(BarcodePair,through=BarcodePairBarcodePlate)
### -------------------------------------------------------------------------------------
class BarcodePairBarcodePlate(models.Model):
    pair = models.ForeignKey(BarcodePair)
    plate = models.ForeignKey(BarcodePlate)
    well = models.CharField(max_length=3)

    class Meta:
        unique_together = (("plate", "well"),)



### -------------------------------------------------------------------------------------
class Panel(models.Model):#collection of targets
    name = models.CharField(max_length=50)
    targets = models.ManyToManyField(TargetEnrichment, related_name='panels')
    #TODO: add comment field
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Primers Multiplex
### -------------------------------------------------------------------------------------
class PrimersMultiplex(models.Model):
    name = models.CharField(max_length=20)
    primers = models.ManyToManyField(TargetEnrichment)
    physical_locations = fields.GenericRelation(SampleLocation,
                               content_type_field='content_type',
                               object_id_field='object_id')

    def __unicode__(self):
        return self.name


