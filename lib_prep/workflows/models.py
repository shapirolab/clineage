from django.contrib.auth.models import User
from django.core.urlresolvers import reverse

from django.db import models
from django.contrib.contenttypes import fields
from linapp.models import Protocol
from wet_storage.models import SampleLocation
from sampling.models import Cell
from primers.parts.models import DNABarcode1, DNABarcode2

class BarcodePair(models.Model):
    left = models.ForeignKey(DNABarcode1)
    right = models.ForeignKey(DNABarcode2)

class CellContentType(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name


class CellContentProtocol(Protocol):
    pass

class CellContent(models.Model):  # aka DNA
    # parent = models.ForeignKey('CellContent', null=True, blank=True)
    cell = models.ForeignKey(Cell)
    # panel = models.ForeignKey(Panel, null=True, blank=True)
    type = models.ForeignKey(CellContentType)
    name = models.CharField(max_length=50, null=True, blank=True)
    protocol = models.ForeignKey(CellContentProtocol, null=True, blank=True)
    # seq_ready = models.BooleanField(default=False)
    user = models.ForeignKey(User, null=True, blank=True)
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


class WorkFlowCell(models.Model): # cell + barcode
    content = models.ForeignKey(CellContent)
    barcodes = models.ForeignKey(BarcodePair)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')
