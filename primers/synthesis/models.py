
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

from primers.parts.models import IlluminaReadingAdaptor1ForTail, \
    IlluminaReadingAdaptor1ForHead, IlluminaReadingAdaptor2ForTail, \
    IlluminaReadingAdaptor2ForHead, IlluminaFlowCellAdaptor2, \
    IlluminaFlowCellAdaptor1, DNABarcode1, DNABarcode2, \
    PadlockAmplificationPlusPrimerPart1, PadlockAmplificationPlusPrimerPart2, \
    PadlockAmplificationMinusPrimerPart1, \
    PadlockAmplificationMinusPrimerPart2, \
    IlluminaReadingAdaptor1, IlluminaReadingAdaptor2, Backbone
from targeted_enrichment.planning.models import UGSPlus, UGSMinus, \
    RestrictionEnzyme
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

    def __str__(self):
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
        return self.iraft.sequence

class PCR1PlusPrimer(PCR1TailMixin,TargetedPlusPrimer):
    """
    Primer that prepends (part of) the Illumina Read Adaptor to a targeted
    amplicon, binding at the UGS leading to our target.
    """
    iraft = models.ForeignKey(IlluminaReadingAdaptor1ForTail)

class PCR1MinusPrimer(PCR1TailMixin,TargetedMinusPrimer):
    """
    Primer that appends (part of) the Illumina Read Adaptor to a targeted
    amplicon, binding at the (rev-comp) UGS following our target.
    """
    iraft = models.ForeignKey(IlluminaReadingAdaptor2ForTail)

class PCR1WithCompanyTagTailMixin(object):
    @property
    def tail(self):
        return self.iraft.sequence + DNA(self.tag)

class PCR1WithCompanyTagPlusPrimer(PCR1WithCompanyTagTailMixin,TargetedPlusPrimer):
    """
    Primer that prepends (part of) the Illumina Read Adaptor, and a random tag
    for identifying the company producing it, to a targeted amplicon, binding
    at the UGS leading to our target.
    """
    tag = models.CharField(max_length=1) # DNAField
    iraft = models.ForeignKey(IlluminaReadingAdaptor1ForTail)

class PCR1WithCompanyTagMinusPrimer(PCR1WithCompanyTagTailMixin,TargetedMinusPrimer):
    """
    Primer that appends (part of) the Illumina Read Adaptor, and a random tag
    for identifying the company producing it, to a targeted amplicon, binding
    at the (rev-comp) UGS following our target.
    """
    tag = models.CharField(max_length=1) # DNAField
    iraft = models.ForeignKey(IlluminaReadingAdaptor2ForTail)

class PCR2Mixin(object):
    @property
    def head(self):
        # NOTE: this is not quite accurate...
        return self.irafh.sequence

    @property
    def tail(self):
        return self.ifca.sequence + self.barcode.sequence

class PCR2PlusPrimer(PCR2Mixin,BasePrimer,PlusStrandMixin):
    """
    Primer that binds at (some of) the Illumina Read Adaptor, prepends the
    rest, the DNA barcode, and the Illumina FlowCell Adaptors.
    """
    irafh = models.ForeignKey(IlluminaReadingAdaptor1ForHead)
    barcode = models.ForeignKey(DNABarcode1)
    ifca = models.ForeignKey(IlluminaFlowCellAdaptor1)

class PCR2MinusPrimer(PCR2Mixin,BasePrimer,MinusStrandMixin):
    """
    Primer that binds at (some of) the Illumina Read Adaptor, appends the
    rest, the DNA barcode, and the Illumina FlowCell Adaptors.
    """
    irafh = models.ForeignKey(IlluminaReadingAdaptor2ForHead)
    barcode = models.ForeignKey(DNABarcode2)
    ifca = models.ForeignKey(IlluminaFlowCellAdaptor2)

class BasePadlockPrep(models.Model):
    """
    A piece of DNA that selectively binds to the template in two locations,
    and by filling in the gap gets a circular DNA containing the targeted
    region.
    Binds to the amplicon at left_region and right_region, and prepends
    left_tail and append right_tail to the new amplicon.
    Contains generic adaptors on each end, for amplification from oligomix.
    """
    name = models.CharField(max_length=50)
    left_amp_primer_part1 = models.ForeignKey(PadlockAmplificationPlusPrimerPart1)
    left_amp_primer_part2 = models.ForeignKey(PadlockAmplificationPlusPrimerPart2)
    right_amp_primer_part1 = models.ForeignKey(PadlockAmplificationMinusPrimerPart1)
    right_amp_primer_part2 = models.ForeignKey(PadlockAmplificationMinusPrimerPart2)
    restriction_enzyme = models.ForeignKey(RestrictionEnzyme)
    # padlock = models.ForeignKey(BasePadlock)
    physical_locations = fields.GenericRelation(SampleLocation,
                                             content_type_field='content_type',
                                             object_id_field='object_id')


    @property
    def sequence(self):
        return self.left_amp_primer_part1.sequence + \
            self.restriction_enzyme.sequence + \
            self.left_amp_primer_part2.sequence + \
            self.padlock.sequence + \
            self.right_amp_primer_part2.sequence + \
            self.restriction_enzyme.sequence.rev_comp() + \
            self.right_amp_primer_part1.sequence

    class Meta:
        abstract = True


class OM6Padlock(models.Model):
    left_ugs = models.ForeignKey(UGSPlus)
    right_ugs = models.ForeignKey(UGSMinus)
    ira1ft = models.ForeignKey(IlluminaReadingAdaptor1ForTail)
    ira2ft = models.ForeignKey(IlluminaReadingAdaptor2ForTail)
    backbone = models.ForeignKey(Backbone)
    umi_length = models.PositiveSmallIntegerField()

    @property
    def sequence(self):
        return self.right_ugs.ref_sequence + \
            DNA.umi(self.umi_length) + \
            self.ira2ft.ref_sequence + \
            self.backbone + \
            self.ira1ft.ref_sequence + \
            DNA.umi(self.umi_length) + \
            self.left_ugs.ref_sequence


class OM6Prep(BasePadlockPrep):
    padlock = models.ForeignKey(OM6Padlock)
