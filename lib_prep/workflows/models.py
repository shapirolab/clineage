from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.db import models
from django.contrib.contenttypes import fields

from model_utils.managers import InheritanceManager

from linapp.models import Protocol
from wet_storage.models import SampleLocation
from sampling.models import Cell
from primers.parts.models import DNABarcode1, DNABarcode2
from lib_prep.multiplexes.models import PCR1MultiplexCollection

class BarcodePair(models.Model):
    left = models.ForeignKey(DNABarcode1)
    right = models.ForeignKey(DNABarcode2)

#class CellContentType(models.Model):
    #name = models.CharField(max_length=50)

    #def __unicode__(self):
        #return self.name


class CellContentProtocol(Protocol):
    pass

class AmplifiedContent(models.Model):  # aka DNA
    # parent = models.ForeignKey('AmplifiedContent', null=True, blank=True)
    cell = models.ForeignKey(Cell)
    # panel = models.ForeignKey(Panel, null=True, blank=True)
    #type = models.ForeignKey(CellContentType)
    name = models.CharField(max_length=50, null=True, blank=True)
    protocol = models.ForeignKey(CellContentProtocol, null=True, blank=True)
    # seq_ready = models.BooleanField(default=False)
    #user = models.ForeignKey(User, null=True, blank=True)
    comment = models.TextField()
    physical_locations = fields.GenericRelation(SampleLocation,
                               content_type_field='content_type',
                               object_id_field='object_id')

    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.cell), self.name)

    def get_absolute_url(self):
        return reverse('cell_content_detail', kwargs={'pk': self.pk})

    def autoname(self):
        if self.seq_ready:
            return "%s_%s_%s" % self.cell.name, self.protocol.initials, 'seqready'
        return "%s_%s" % self.cell.name, self.protocol.initials

    # def experiment(self):
    #    print self.cell.experiment.values('id').annotate(experiment_count=Count('id')).order_by('-experiment_count')
    #    return Experiment.objects.get(id = self.cell.experiment.values('id').annotate(experiment_count=Count('id')).order_by('-experiment_count')[0]['id'])

class Library(models.Model):
    name = models.CharField(max_length=50)

    @property
    def barcoded_contents(self):
        raise NotImplementedError()

    @property
    def unwrappers(self):
        """
        Return an iterator (favorably, a QuerySet) for unwrappers which might
        match this cell.
        """
        raise NotImplementedError()

    objects = InheritanceManager()

class BarcodedContent(models.Model): # cell + barcode
    barcodes = models.ForeignKey(BarcodePair)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')

    @property
    def cell(self):
        raise NotImplementedError()

class MagicalPCR1Library(Library):
    mpx_collection = models.ForeignKey(PCR1MultiplexCollection)
    
    @property
    def barcoded_contents(self):
        return self.magicalpcr1barcodedcontent_set.all()

    @property
    def unwrappers(self):
        #TODO: make nice and queryful.
        for mpx in self.mpx_collection.mpxs.all():
            for primer in mpx.primers.all():
                yield primer.pcr1unwrapper

class MagicalPCR1BarcodedContent(BarcodedContent):
    content = models.ForeignKey(AmplifiedContent)
    library = models.ForeignKey(MagicalPCR1Library)

    @property
    def cell(self):
        return self.content.cell
