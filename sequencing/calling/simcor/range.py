from django.db import models
from sequencing.calling.simcor.simulation_spaces import get_far_apart_highest_peaks
from sequencing.calling.range import MultiAlleleMixin, AllelesRangeMixin, ProportionalAllelesMixin, \
    ProportionsRangeMixin, BoundProportionsRangeMixin
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


class AllelesCyclesRangeMixin(AllelesRangeMixin, CyclesRangeMixin, BaseAllelesCyclesRange):
    """
    Brute forcing over all monoallelic options at all cycles
    """
    pass


class FullRangeBiMixin(MultiAlleleMixin, CyclesRangeMixin, BaseAllelesCyclesRange):
    """
    Brute forcing over all bi-allelic options at all cycles
    """
    pass


class ProportionalAllelesCyclesRangeMixin(ProportionsRangeMixin, ProportionalAllelesMixin,
                                          CyclesRangeMixin, BaseAllelesCyclesRange):
    """
    Brute forcing over all biallelic options at all proportions and all cycles
    """
    pass


class BoundProportionalAllelesCyclesRangeMixin(BoundProportionsRangeMixin, ProportionalAllelesMixin,
                                               CyclesRangeMixin, BaseAllelesCyclesRange):
    """
    Brute forcing over all biallelic options at all proportions and all cycles
    """
    pass


class HighestPeaksMixin(ProportionalAllelesMixin):

    @classmethod
    def highest_peaks(cls, hist):
        return get_far_apart_highest_peaks(
            hist=hist,
            allele_number=cls.allele_number,
            minimal_distance_between_peaks=1,
        )


class AllelesRangeFromPointModelMixin(models.Model):
    range_from_point = models.PositiveSmallIntegerField()

    class Meta:
        abstract = True


class HighestPeaksModelMixin(AllelesRangeFromPointModelMixin, HighestPeaksMixin):

    def alleles_by_hist(self, hist):
        points = self.highest_peaks(hist)
        yield from itertools.chain(
            range(point-self.range_from_point, point+self.range_from_point+1)
            for point in points
        )


class PeaksMinimalDistanceModelMixin(models.Model, HighestPeaksMixin):
    minimal_seeds_distance = models.PositiveSmallIntegerField()

    class Meta:
        abstract = True

    def highest_peak(self, hist):
        return get_far_apart_highest_peaks(
            hist=hist,
            allele_number=self.allele_number,
            minimal_distance_between_peaks=self.minimal_seeds_distance,
        )


class HighestPeaksRangeMixin(PeaksMinimalDistanceModelMixin, HighestPeaksModelMixin):

    pass
