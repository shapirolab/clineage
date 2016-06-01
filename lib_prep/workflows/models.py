from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.db import models
from django.contrib.contenttypes import fields

from model_utils.managers import InheritanceManager

from linapp.models import Protocol
from wet_storage.models import SampleLocation
from sampling.models import Cell
from primers.parts.models import DNABarcode1, DNABarcode2
from primers.synthesis.models import PCR2PlusPrimer, PCR2MinusPrimer
from lib_prep.multiplexes.models import PCR1Panel, OM6Panel
from targeted_enrichment.reagents.models import TwoPrimersStrMixin

class BarcodePair(models.Model):
    left = models.ForeignKey(DNABarcode1)
    right = models.ForeignKey(DNABarcode2)

    def __str__(self):
        return "{}, {}".format(self.left, self.right)


class PCR2PrimerPairReagent(TwoPrimersStrMixin, models.Model):
    bp = models.ForeignKey(BarcodePair)  # TODO: maybe kill?
    left_primer = models.ForeignKey(PCR2PlusPrimer)
    right_primer = models.ForeignKey(PCR2MinusPrimer)
    physical_locations = fields.GenericRelation(SampleLocation,
                                                 content_type_field='content_type',
                                                 object_id_field='object_id')


class CellContentProtocol(Protocol):
    pass


class AmplifiedContent(models.Model):  # aka DNA
    cell = models.ForeignKey(Cell)
    name = models.CharField(max_length=50, null=True, blank=True)
    protocol = models.ForeignKey(CellContentProtocol, null=True, blank=True)
    comment = models.TextField()
    physical_locations = fields.GenericRelation(SampleLocation,
                               content_type_field='content_type',
                               object_id_field='object_id')

    def __str__(self):
        return '{}>{}'.format(self.cell, self.name)

    def get_absolute_url(self):
        return reverse('cell_content_detail', kwargs={'pk': self.pk})


class Library(models.Model):
    name = models.CharField(max_length=50)

    # FIXME
    @property
    def subclass(self):
        return Library.objects.get_subclass(id=self.id)

    @property
    def barcoded_contents(self):
        raise NotImplementedError()

    @property
    def amplicons(self):
        """
        Return an iterator (favorably, a QuerySet) for amplicons which might
        match this cell.
        """
        raise NotImplementedError()

    objects = InheritanceManager()

    def __str__(self):
        return self.name


class BarcodedContent(models.Model): # cell + barcode
    barcodes = models.ForeignKey(BarcodePair)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')

    objects = InheritanceManager()

    # FIXME
    @property
    def subclass(self):
        return BarcodedContent.objects.get_subclass(id=self.id)

    @property
    def amplified_content(self):
        raise NotImplementedError()

    def __str__(self):
        return '{}/({})'.format(self.amplified_content, self.barcodes)


class UnsupportedLibrary(Library):

    @property
    def barcoded_contents(self):
        return self.unsupportedbarcodedcontent_set.all()

    @property
    def amplicons(self):
        return ()


class UnsupportedBarcodedContent(BarcodedContent):
    content = models.ForeignKey(AmplifiedContent)
    library = models.ForeignKey(UnsupportedLibrary)

    @property
    def amplified_content(self):
        return self.content


class MagicalPCR1Library(Library):
    panel = models.ForeignKey(PCR1Panel)
    # magicalpcr1barcodedcontent_set is a related field

    @property
    def barcoded_contents(self):
        return self.magicalpcr1barcodedcontent_set.all()

    @property
    def amplicons(self):
        #TODO: make nice and queryful.
        for mpx in self.panel.mpxs.all():
            for ter in mpx.ters.select_subclasses():
                yield ter.amplicon


class MagicalPCR1BarcodedContent(BarcodedContent):
    content = models.ForeignKey(AmplifiedContent)
    library = models.ForeignKey(MagicalPCR1Library)

    @property
    def amplified_content(self):
        return self.content


class MagicalOM6Library(Library):
    panel = models.ForeignKey(OM6Panel)
    # magicalom6barcodedcontent_set is a related field

    @property
    def barcoded_contents(self):
        return self.magicalom6barcodedcontent_set.all()

    @property
    def amplicons(self):
        #TODO: make nice and queryful.
        for mix in self.panel.mixs.all():
            for ter in mix.ters.select_subclasses():
                yield ter.amplicon


class MagicalOM6BarcodedContent(BarcodedContent):
    content = models.ForeignKey(AmplifiedContent)
    library = models.ForeignKey(MagicalOM6Library)

    @property
    def amplified_content(self):
        return self.content
