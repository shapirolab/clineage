from django.db import models
from sequencing.calling.simcor.simulation_spaces import get_far_apart_highest_peaks
from sequencing.calling.range import MultiAlleleMixin, AllelesRangeMixin, ProportionalAllelesMixin, \
    BoundProportionalAllelesMixin
from sequencing.calling.simcor.models_common import CyclesRangeMixin
import itertools


class BaseAllelesCyclesRange(object):

    @property
    def alleles(self):
        raise NotImplemented

    @property
    def cycles(self):
        raise NotImplemented

    @property
    def alleles_and_cycles(self):
        yield from itertools.product(self.alleles, self.cycles)


class AllelesCyclesRangeMixin(BaseAllelesCyclesRange, AllelesRangeMixin, CyclesRangeMixin):
    """
    Brute forcing over all monoallelic options at all cycles
    """
    pass


class FullRangeBiMixin(BaseAllelesCyclesRange, MultiAlleleMixin, CyclesRangeMixin):
    """
    Brute forcing over all bi-allelic options at all cycles
    """
    pass


class ProportionalAllelesCyclesRangeMixin(BaseAllelesCyclesRange, ProportionalAllelesMixin, CyclesRangeMixin):
    """
    Brute forcing over all biallelic options at all proportions and all cycles
    """
    pass


class BoundProportionalAllelesCyclesRangeMixin(BaseAllelesCyclesRange, BoundProportionalAllelesMixin, CyclesRangeMixin):
    """
    Brute forcing over all biallelic options at all proportions and all cycles
    """
    pass


class HighestPeakMixin(CyclesRangeMixin, MultiAlleleMixin):

    @classmethod
    def highest_peak(cls, hist):
        return get_far_apart_highest_peaks(
            hist=hist,
            allele_number=cls.allele_number,
            minimal_distance_between_peaks=1,
        )


class PeaksMinimalDistanceModelMixin(models.Model, HighestPeakMixin):
    minimal_seeds_distance = models.PositiveSmallIntegerField()

    class Meta:
        abstract = True

    def highest_peak(self, hist):
        return get_far_apart_highest_peaks(
            hist=hist,
            allele_number=self.allele_number,
            minimal_distance_between_peaks=self.minimal_seeds_distance,
        )
