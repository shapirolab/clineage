
from django.db import models
from django.contrib.contenttypes import fields

# Multiplexing Read 1 Sequencing Primer
# 5'
# ACACTCTTTCCCTACACGACGCTCTTCCGATCT
# Multiplexing Index Read Sequencing Primer
# 5'
# GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# Multiplexing Read 2 Sequencing Primer
# 5'
# GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
#
#
# AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT[INSERT][A]GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG

from primers.parts.models import IlluminaReadingAdaptor1Cuts, \
    IlluminaReadingAdaptor2Cuts, IlluminaFlowCellAdaptor2, \
    IlluminaFlowCellAdaptor1, DNABarcode1, DNABarcode2, \
    PadlockAmplificationPlusPrimer, PadlockAmplificationMinusPrimer
from targeted_enrichment.planning.models import UGSPlus, UGSMinus
from primers.strand import BaseStrandMixin, MinusStrandMixin, PlusStrandMixin
from misc.dna import DNA

from wet_storage.models import SampleLocation

class BasePrimer(models.Model,BaseStrandMixin):
    """
    Any DNA that can act as a primer for PCR.
    Binds to the amplicon at the head (amplicon should begin with head)
    and optionally prepends the tail to the new amplicon.
    """
    name = models.CharField(max_length=50)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')

    @property
    def head(self):
        raise NotImplementedError()

    @property
    def tail(self):
        raise NotImplementedError()

    @property
    def sequence(self):
        return self.tail + self.head

    class Meta:
        abstract = True

    def __unicode__(self):
        return self.name

class TargetedHeadMixin(object):
    @property
    def head(self):
        return self.ugs.sequence

class TargetedPlusPrimer(TargetedHeadMixin,BasePrimer,PlusStrandMixin):
    ugs = models.ForeignKey(UGSPlus)

    class Meta:
        abstract = True

class TargetedMinusPrimer(TargetedHeadMixin,BasePrimer,MinusStrandMixin):
    ugs = models.ForeignKey(UGSMinus)

    class Meta:
        abstract = True

class NoTailMixin(object):
    @property
    def tail(self):
        return DNA("")

class TargetedNoTailPlusPrimer(NoTailMixin,TargetedPlusPrimer):
    pass

class TargetedNoTailMinusPrimer(NoTailMixin,TargetedMinusPrimer):
    pass

class PCR1TailMixin(object):
    @property
    def tail(self):
        return self.irac.primer1tail

class PCR1PlusPrimer(PCR1TailMixin,TargetedPlusPrimer):
    """
    Primer that prepends (part of) the Illumina Read Adaptor to a targeted
    amplicon, binding at the UGS leading to our target.
    """
    irac = models.ForeignKey(IlluminaReadingAdaptor1Cuts)

class PCR1MinusPrimer(PCR1TailMixin,TargetedMinusPrimer):
    """
    Primer that appends (part of) the Illumina Read Adaptor to a targeted
    amplicon, binding at the (rev-comp) UGS following our target.
    """
    irac = models.ForeignKey(IlluminaReadingAdaptor2Cuts)

class PCR1WithCompanyTagTailMixin(object):
    @property
    def tail(self):
        return self.irac.primer1tail + DNA(self.tag)

class PCR1WithCompanyTagPlusPrimer(PCR1WithCompanyTagTailMixin,TargetedPlusPrimer):
    """
    Primer that prepends (part of) the Illumina Read Adaptor, and a random tag
    for identifying the company producing it, to a targeted amplicon, binding
    at the UGS leading to our target.
    """
    tag = models.CharField(max_length=1) # DNAField
    irac = models.ForeignKey(IlluminaReadingAdaptor1Cuts)

class PCR1WithCompanyTagMinusPrimer(PCR1WithCompanyTagTailMixin,TargetedMinusPrimer):
    """
    Primer that appends (part of) the Illumina Read Adaptor, and a random tag
    for identifying the company producing it, to a targeted amplicon, binding
    at the (rev-comp) UGS following our target.
    """
    tag = models.CharField(max_length=1) # DNAField
    irac = models.ForeignKey(IlluminaReadingAdaptor2Cuts)

class PCR2Mixin(object):
    @property
    def head(self):
        return self.irac.overlap

    @property
    def tail(self):
        return self.ifca.sequence + self.barcode.sequence + self.irac.primer2tail

class PCR2PlusPrimer(PCR2Mixin,BasePrimer,PlusStrandMixin):
    """
    Primer that binds at (some of) the Illumina Read Adaptor, prepends the
    rest, the DNA barcode, and the Illumina FlowCell Adaptors.
    """
    irac = models.ForeignKey(IlluminaReadingAdaptor1Cuts)
    barcode = models.ForeignKey(DNABarcode1)
    ifca = models.ForeignKey(IlluminaFlowCellAdaptor1)

class PCR2MinusPrimer(PCR2Mixin,BasePrimer,MinusStrandMixin):
    """
    Primer that binds at (some of) the Illumina Read Adaptor, appends the
    rest, the DNA barcode, and the Illumina FlowCell Adaptors.
    """
    irac = models.ForeignKey(IlluminaReadingAdaptor2Cuts)
    barcode = models.ForeignKey(DNABarcode2)
    ifca = models.ForeignKey(IlluminaFlowCellAdaptor2)

class BaseShortPadlock(models.Model):
    """
    A piece of DNA that selectively binds to the template in two locations,
    and by filling in the gap gets a circular DNA containing the targeted
    region.
    Binds to the amplicon at left_region and right_region, and prepends
    left_tail and append right_tail to the new amplicon.
    Contains generic adaptors on each end, for amplification from oligomix.
    """
    name = models.CharField(max_length=50)
    left_amp_primer = models.ForeignKey(PadlockAmplificationPlusPrimer)
    right_amp_primer = models.ForeignKey(PadlockAmplificationMinusPrimer)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')

    @property
    def left_region(self):
        raise NotImplementedError()

    @property
    def right_region(self):
        raise NotImplementedError()

    @property
    def left_tail(self):
        raise NotImplementedError()

    @property
    def right_tail(self):
        raise NotImplementedError()

    @property
    def backbone(self):
        raise NotImplementedError()

    @property
    def sequence(self):
        return self.left_amp_primer.sequence + \
            self.right_region + self.right_tail + \
            self.backbone + \
            self.left_tail + self.left_region + \
            self.right_amp_primer.sequence

    class Meta:
        abstract = True

class BaseLongPadlockTemplate(models.Model):
    """
    The original template from which a long padlock is generated, using stock
    DNA as the backbone.
    Contains the left and right targeting regions separated by a restriction
    site, which is digested after annealing with the backbone.
    """
    name = models.CharField(max_length=50)
    left_amp_primer = models.ForeignKey(PadlockAmplificationPlusPrimer)
    right_amp_primer = models.ForeignKey(PadlockAmplificationMinusPrimer)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')

    @property
    def left_region(self):
        raise NotImplementedError()

    @property
    def right_region(self):
        raise NotImplementedError()

    @property
    def left_tail(self):
        raise NotImplementedError()

    @property
    def right_tail(self):
        raise NotImplementedError()

    @property
    def filler(self):
        raise NotImplementedError()

    @property
    def left_backbone_adaptor(self):
        raise NotImplementedError()

    @property
    def right_backbone_adaptor(self):
        raise NotImplementedError()

    @property
    def sequence(self):
        return self.right_backbone_adaptor + \
            self.left_tail + self.left_region + \
            self.right_amp_primer.sequence + \
            self.filler + \
            self.left_amp_primer.sequence + \
            self.right_region + self.right_tail + \
            self.left_backbone_adaptor

    class Meta:
        abstract = True
