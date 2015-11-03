from django.db import models
from django.contrib.contenttypes import fields
from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.conf import settings

from misc.models import Taxa
from wet_storage.models import SampleLocation


### -------------------------------------------------------------------------------------
### Types and descriptors
### -------------------------------------------------------------------------------------
class GeneticBackground(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Organ(models.Model):
    name = models.CharField(max_length=50)
    #rank = models.IntegerField()
    #parent = models.IntegerField()

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Tissue(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class SampleComposition(models.Model):#e.g. single cell or bulk
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------

class SampleStatus(models.Model):
    name = models.CharField(max_length=50)

    class Meta:
        verbose_name = 'Sample status'
        verbose_name_plural = 'Samples status'

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------

### -------------------------------------------------------------------------------------
class Coordinates(models.Model):#XYZ coordinates of laser capture.
    x = models.DecimalField(max_digits=10, decimal_places=4)
    y = models.DecimalField(max_digits=10, decimal_places=4)
    z = models.DecimalField(max_digits=10, decimal_places=4)
    #TODO: discuss Z vs slides
    def __unicode__(self):
        return '(x={0},y={1},z={2})'.format(self.x,self.y,self.z)
### -------------------------------------------------------------------------------------
class FACSMarker(models.Model):
    name = models.CharField(max_length=50)
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Location(models.Model):  # Freetext location
    name = models.CharField(max_length=50)

    def __unicode__(self):#TODO: exapnd
        return self.name
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Sampling Hierarchy
### -------------------------------------------------------------------------------------
class Individual(models.Model):
    GENDER = (('M', 'Male'), ('F', 'Female'),)
    taxa = models.ForeignKey(Taxa)
    name = models.CharField(max_length=50)
    sex = models.CharField(max_length=1, choices=GENDER)
    born = models.DateTimeField(null=True, blank=True)
    comment = models.TextField(null=True, blank=True)
    background = models.ForeignKey(GeneticBackground, null=True, blank=True)
    location = models.ForeignKey(Location, null=True, blank=True)
    sacrificed = models.DateTimeField(null=True, blank=True)
    partner = models.ForeignKey(User, null=True, blank=True)

    def __unicode__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('individual_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class ExtractionEvent(models.Model):
    individual = models.ForeignKey(Individual)
    name = models.CharField(max_length=100)
    comment = models.TextField(null=True, blank=True)
    date = models.DateTimeField()
    location = models.ForeignKey(Location, null=True, blank=True)
    user_performed = models.ForeignKey(User, related_name='+')
    user_documented = models.ForeignKey(User, related_name='+')

    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.individual), self.name)

    def get_absolute_url(self):
        return reverse('extraction_event_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class Extraction(models.Model):
    extraction_event = models.ForeignKey(ExtractionEvent)
    name = models.CharField(max_length=50)
    date = models.DateTimeField(null=True, blank=True)
    organ = models.ForeignKey(Organ)
    tissue = models.ForeignKey(Tissue)
    comment = models.TextField(null=True, blank=True)
    physical_locations = fields.GenericRelation('SampleLocation',
                               content_type_field='content_type',
                               object_id_field='object_id')
    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.extraction_event), self.name)

    def get_absolute_url(self):
        return reverse('extraction_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
def sampling_event_path(self, filename):
    return settings.MEDIA_ROOT+'/CLineageFiles/sampling_events/'+self.name+'/'+filename
### -------------------------------------------------------------------------------------
class SamplingEvent(models.Model):
    name = models.CharField(max_length=50)
    extraction = models.ForeignKey(Extraction)
    date = models.DateField(null=True, blank=True)
    user = models.ForeignKey(User, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)
    attachment = models.FileField(upload_to=sampling_event_path, null=True, blank=True)
    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.extraction), self.name)

    def get_absolute_url(self):
        return reverse('samplingevent_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class FACS(SamplingEvent):
    marker = models.ForeignKey(FACSMarker)

    def get_absolute_url(self):
        return reverse('facs_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class LaserCapture(SamplingEvent):
    coordinates = models.ForeignKey(Coordinates, null=True, blank=True)  # every cell has distinct coordinates.

    def get_absolute_url(self):
        return reverse('lasercapture_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class CellSelector(SamplingEvent):
    coordinates = models.ForeignKey(Coordinates, null=True, blank=True)  # every cell has distinct coordinates.

    def get_absolute_url(self):
        return reverse('cellselector_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class Cell(models.Model):
    individual = models.ForeignKey(Individual)
    sampling = models.ForeignKey(SamplingEvent, null=True)
    name = models.CharField(max_length=50)
    #experiment = models.ManyToManyField(Experiment, related_name='cells', null=True, blank=True)
    composition = models.ForeignKey(SampleComposition)#single cell or bulk
    status = models.ForeignKey(SampleStatus, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)
    classification = models.CharField(max_length=50, null=True, blank=True)
    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.sampling), self.name)

    def get_absolute_url(self):
        return reverse('cell_detail', kwargs={'pk': self.pk})

### -------------------------------------------------------------------------------------
