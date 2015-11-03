
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
    IlluminaFlowCellAdaptor1, DNABarcode1, DNABarcode2
from targeted_enrichment.planning.models import UGSPlus, UGSMinus
from primers.strand import BaseStrandMixin, MinusStrandMixin, PlusStrandMixin

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
        return ""

    @property
    def tail(self):
        return ""

    @property
    def sequence(self):
        return self.tail + self.head

    class Meta:
        abstract = True

class TargetedHeadMixin(object):
    @property
    def head(self):
        return self.ugs.sequence

class TargetedNoTailPlusPrimer(TargetedHeadMixin,BasePrimer,PlusStrandMixin):
    ugs = models.ForeignKey(UGSPlus)

class TargetedNoTailMinusPrimer(TargetedHeadMixin,BasePrimer,MinusStrandMixin):
    ugs = models.ForeignKey(UGSMinus)

class PCR1TailMixin(object):
    @property
    def tail(self):
        return self.irac.primer1tail

class PCR1PlusPrimer(TargetedHeadMixin,PCR1TailMixin,BasePrimer,PlusStrandMixin):
    """
    Primer that prepends (part of) the Illumina Read Adaptor to a targeted
    amplicon, binding at the UGS leading to our target.
    """
    ugs = models.ForeignKey(UGSPlus)
    irac = models.ForeignKey(IlluminaReadingAdaptor1Cuts)

class PCR1MinusPrimer(TargetedHeadMixin,PCR1TailMixin,BasePrimer,MinusStrandMixin):
    """
    Primer that appends (part of) the Illumina Read Adaptor to a targeted
    amplicon, binding at the (rev-comp) UGS following our target.
    """
    ugs = models.ForeignKey(UGSMinus)
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

