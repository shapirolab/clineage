from django.db import models
from sequencing.calling.simcor.hist_analysis import get_far_apart_highest_peaks
from sequencing.calling.range import MultiAlleleMixin, AllelesRangeMixin, ProportionalAllelesMixin, \
    ProportionsRangeMixin, BoundProportionsRangeMixin, BaseAlleleNumber
from sequencing.calling.simcor.models_common import CyclesRangeMixin
import itertools
import functools
import decimal


class BaseAllelesCyclesRangeMixin(object):
    """
    Base Mixin for keeping alleles and cycles for current calling schema
    """

    @property
    def alleles(self):
        raise NotImplementedError

    @property
    def cycles(self):
        raise NotImplementedError

    @property
    def alleles_and_cycles(self):
        yield from itertools.product(self.alleles, self.cycles)


class AllelesCyclesRangeMixin(AllelesRangeMixin, CyclesRangeMixin, BaseAllelesCyclesRangeMixin):
    """
    Brute forcing over all monoallelic options at all cycles
    """
    pass


class FullRangeBiMixin(MultiAlleleMixin, CyclesRangeMixin, BaseAllelesCyclesRangeMixin):
    """
    Brute forcing over all bi-allelic options at all cycles
    """
    pass


class ProportionalAllelesCyclesRangeMixin(ProportionsRangeMixin, ProportionalAllelesMixin,
                                          CyclesRangeMixin, BaseAllelesCyclesRangeMixin):
    """
    Brute forcing over all biallelic options at all proportions and all cycles
    """
    pass


class BoundProportionalAllelesCyclesRangeMixin(BoundProportionsRangeMixin, ProportionalAllelesMixin,
                                               CyclesRangeMixin, BaseAllelesCyclesRangeMixin):
    """
    Brute forcing over all biallelic options at all proportions and all cycles
    """
    pass


def pairwise(iterable):
    """
    itertools recipe
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


@functools.lru_cache(maxsize=4096)  # apparantly this can be heavy
def get_proportion_bounds(diff, longer_ms_len, cycle=20, length_sensitivity=decimal.Decimal(0.21),
                          diff_sensetivity=decimal.Decimal(0.65), eps=decimal.Decimal(0.1)):
    """
    This function calculates the valid proportion bounds according to the difference between the alleles and length of the longer allele
    *Defaults are based on old AC experiment
    """
    equilibrium = decimal.Decimal(0.5)
    distance_from_equilibrium = min(
        ((diff ** diff_sensetivity / (length_sensitivity * longer_ms_len)) ** 3),
        equilibrium)
    if distance_from_equilibrium < eps:
        return equilibrium, equilibrium
    return equilibrium-distance_from_equilibrium, equilibrium+distance_from_equilibrium


@functools.lru_cache(maxsize=16384)
def contains_excluded_proportions(alleles_and_proportions, cycle=20, length_sensitivity=decimal.Decimal(0.21), diff_sensetivity=decimal.Decimal(0.65)):
    """
    This function checks for each allele and proportion whether it's in calculated bounds
    *Defaults are based on old AC experiment
    """
    mono_cases = [decimal.Decimal('1.0'), decimal.Decimal('0.0')]
    for (a1, p1), (a2, p2) in pairwise(alleles_and_proportions):
        lower_bound, upper_bound = get_proportion_bounds(abs(a1-a2), max(a1, a2), cycle=cycle, length_sensitivity=length_sensitivity, diff_sensetivity=diff_sensetivity)
        if not all(lower_bound < p < upper_bound or p in mono_cases for p in [p1, p2]):
            # if any of the proportions (p) is out of bounds
            return True
    return False


def contains_excluded_proportions_wrapper(alleles_and_cycles, length_sensitivity=0.21, diff_sensetivity=0.65):
    alleles_and_proportions, cycle = alleles_and_cycles
    return contains_excluded_proportions(
        alleles_and_proportions,
        # cycle=cycle,  # cycle is not yet implemented in the exclusion function. we drop it to improve caching
        length_sensitivity=length_sensitivity,
        diff_sensetivity=diff_sensetivity)


class ProximityRatioFilteredAllelesCyclesRangeMixin(BaseAllelesCyclesRangeMixin):
    """
    Mixin that sieves alleles_and_cycles using the function contains_excluded_proportions_wrapper.
    These proportions will be used later on as bi-allelic options.
    """

    def prf_filtered(self, alleles_and_cycles_generator):
        yield from itertools.filterfalse(
            functools.partial(
                contains_excluded_proportions_wrapper,
                length_sensitivity=self.length_sensitivity,
                diff_sensetivity=self.diff_sensetivity
            ),
            alleles_and_cycles_generator,
        )

    @property
    def alleles_and_cycles(self):
        yield from self.prf_filtered(
            super().alleles_and_cycles,
        )


class ProximityRatioFilteredProportionalAllelesCyclesRangeMixin(ProportionsRangeMixin, ProportionalAllelesMixin,
                                          CyclesRangeMixin, ProximityRatioFilteredAllelesCyclesRangeMixin):
    """
    Brute forcing over all biallelic options at proportions and cycles bound according to the exclusion function
    """
    pass


class ProximityRatioFilteredBoundProportionalAllelesCyclesRangeMixin(BoundProportionsRangeMixin, ProportionalAllelesMixin,
                                               CyclesRangeMixin, ProximityRatioFilteredAllelesCyclesRangeMixin):
    """
    Brute forcing over all biallelic options at proportions and cycles bound by both hard bounds and
    according to the exclusion function
    """
    pass


class HighestPeaksMixin(BaseAlleleNumber):

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


class PeaksMinimalDistanceModelMixin(models.Model, HighestPeaksMixin):
    minimal_seeds_distance = models.PositiveSmallIntegerField()

    class Meta:
        abstract = True

    def highest_peaks(self, hist):
        return get_far_apart_highest_peaks(
            hist=hist,
            allele_number=self.allele_number,
            minimal_distance_between_peaks=self.minimal_seeds_distance,
        )


class HighestPeaksModelMixin(AllelesRangeFromPointModelMixin, HighestPeaksMixin):

    class Meta:
        abstract = True

    def alleles_by_hist(self, hist):
        points = self.highest_peaks(hist)
        yield from itertools.chain(*[
            range(point-self.range_from_point, point+self.range_from_point+1)
            for point in points]
        )


class HighestPeaksRangeModelMixin(PeaksMinimalDistanceModelMixin, HighestPeaksModelMixin):

    class Meta:
        abstract = True
