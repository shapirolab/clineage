
from django.db import models
from django.contrib.contenttypes import fields
from django.contrib.contenttypes.models import ContentType
from django.core.urlresolvers import reverse

from utils.wells import index2str, str2index


### -------------------------------------------------------------------------------------
### Optional storage mapping
### -------------------------------------------------------------------------------------
class StorageType(models.Model):
    name = models.CharField(max_length=100, blank=True)
    temperature = models.DecimalField(null=True, blank=True, max_digits=5, decimal_places=1)
    def __str__(self):
        return "%s (temp:%s\xB0c)" % (self.name, str(self.temperature))
### -------------------------------------------------------------------------------------
class StorageBox(models.Model):
    code = models.AutoField(primary_key=True)
    storage_type = models.ForeignKey(StorageType)
    name = models.CharField(max_length=100, blank=True)
    barcode = models.CharField(max_length=20, blank=True)
    def __str__(self):
        return "%s (type:%s)" % (self.name, self.storage_type.name)
### -------------------------------------------------------------------------------------
class PlateContext(models.Model): #The plate's context in use. e.g. pcr
    description = models.CharField(max_length=30, blank=True)
    def __str__(self):
        return "%s" %self.description
### -------------------------------------------------------------------------------------
class PlatePlastica(models.Model): #The plate's physical form. e.g. deepwell square
    code = models.AutoField(primary_key=True)
    description = models.CharField(max_length=30, blank=True)
    rows = models.IntegerField(default=8)
    columns = models.IntegerField(default=12)

    @property
    def wells(self):
        return self.rows*self.columns

    def __str__(self):
        return "%s" %self.description
### -------------------------------------------------------------------------------------
class PlateType(models.Model):  # TODO: kill and modify as follows:
    code = models.AutoField(primary_key=True)
    friendly = models.CharField(max_length=100)
    context = models.ForeignKey(PlateContext, null=True)  # TODO: kill and modify to plate inheritance
    plastic = models.ForeignKey(PlatePlastica, null=True)  # TODO: relocate as field of Plate
    def __str__(self):
        return "%s" % self.friendly
### -------------------------------------------------------------------------------------
class PlateFullException(Exception):
    """
    Attempt to position reagent when plate is already full
    """
    pass

class Plate(models.Model):
    code = models.AutoField(primary_key=True)
    type = models.ForeignKey(PlateType)
    name = models.CharField(max_length=200, blank=True)
    barcode = models.CharField(max_length=20, blank=True)
    timestamp = models.DateField(null=True, blank=True)
    state = models.CharField(max_length=20, blank=True)
    lastusedwell = models.CharField(max_length=4, default='A1')
    def __str__(self):
        return self.name
    def get_absolute_url(self):
        return reverse('plate_detail', kwargs={'pk': self.pk})

    def skip_to_first_free_column(self):  # Can't really remember what this function is for. CADMAD heritage.
        first_free_well_index = str2index(self.lastusedwell)
        if first_free_well_index % 8 != 1:
            self.lastusedwell = index2str(first_free_well_index + 8 - (first_free_well_index-1) % 8)

    def place_reagent(self, reagent, reagent_volume=-1, reagent_concentration=-1):
        print(str2index(self.lastusedwell), self.type.plastic.wells)
        if str2index(self.lastusedwell) > self.type.plastic.wells:
            raise PlateFullException
        location = SampleLocation.objects.create(plate=self,
                                                 well=self.lastusedwell,
                                                 reagent=reagent,
                                                 volume=reagent_volume,
                                                 concentration=reagent_concentration)
        self.lastusedwell = index2str(str2index(self.lastusedwell)+1)  # increment first free well
        return location
    place_reagent.PlateFullException = PlateFullException
### -------------------------------------------------------------------------------------
class PlateStorage(models.Model):
    storageBox = models.ForeignKey(StorageBox)
    plate = models.ForeignKey(Plate)
    inner_location = models.CharField(max_length=100, blank=True)
    notes = models.CharField(max_length=250, blank=True)
    def __str__(self):
        return '%s in %s' % (self.plate.name, self.storageBox.name)
### -------------------------------------------------------------------------------------
class SampleLocation(models.Model):
    plate = models.ForeignKey(Plate)
    well = models.CharField(max_length=3, blank=True, db_index=True)
    content_type = models.ForeignKey(ContentType)
    object_id = models.PositiveIntegerField(db_index=True)
    reagent = fields.GenericForeignKey('content_type', 'object_id')
    volume = models.DecimalField(null=True, max_digits=10, decimal_places=3, blank=True)
    concentration = models.DecimalField(null=True, max_digits=10, decimal_places=5, blank=True)
    timestamp = models.DateTimeField(auto_now=True)
    def __str__(self):
        return 'plate %s at %s' % (self.plate.name, self.well)

    class Meta:
        index_together = [
            ["content_type", "object_id"],
        ]
        #unique_together = (
        #    ("plate", "well"),
        #)
