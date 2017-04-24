from django.db import models
from sequencing.calling.models import CallingScheme
from sequencing.calling.simcor.simulation_spaces import mono_sim_hists_space_generator, bi_sim_hists_space_generator,\
    proportional_bi_sim_hists_space_generator
from sequencing.calling.hist_dist import pop_dist_corr_numpy
from sequencing.calling.simcor.calling import call_microsatellite_histogram, get_closest
from sequencing.calling.simcor.models_common import CyclesModelMixin, SimulationsByCycles, MSLengthBoundsModelMixin, \
    ProportionsBoundsModelMixin, ProportionStepModelMixin, BestCorrelationCalledAlleles, \
    BestCorrelationProportionalCalledAlleles, \
    BestCorrelationProportionalHighestPeakCalledAlleles
from sequencing.calling.simcor.range import AllelesCyclesRangeMixin, FullRangeBiMixin, \
    ProportionalAllelesCyclesRangeMixin, BoundProportionalAllelesCyclesRangeMixin, \
    HighestPeaksRangeModelMixin
from itertools import filterfalse


class BestCorrelationCalledAlleleMixin(object):
    @property
    def called_allele_class(self):
        return BestCorrelationCalledAlleles


class BestCorrelationProportionalCalledAlleleMixin(object):
    @property
    def called_allele_class(self):
        return BestCorrelationProportionalCalledAlleles


class BestCorrelationProportionalHighestPeakCalledAlleleMIxin(object):
    @property
    def called_allele_class(self):
        return BestCorrelationProportionalHighestPeakCalledAlleles


class DynamicFilteredHistSpaceMixin(object):
    """
    Filters sim_hists_space by a live histogram object
    """

    def filtered_sim_hists_space(self, hist):
        raise NotImplemented

    def find_best_in_space(self, hist):
        return get_closest(hist, self.filtered_sim_hists_space(hist), self.distance_metric)


class FilterByHistMixin(DynamicFilteredHistSpaceMixin):

    def filtered_sim_hists_space(self, hist):
        """cuts the simulations based on the hist"""
        def allele_in_hist_space(sim_hist):
            alleles_by_hist = self.alleles_by_hist(hist)
            for allele in sim_hist.allele_frozenset:
                if allele not in alleles_by_hist:
                    return True
            return False
        yield from filterfalse(allele_in_hist_space, self.sim_hists_space)


class BaseSimCallingScheme(BestCorrelationCalledAlleleMixin, CallingScheme, MSLengthBoundsModelMixin, CyclesModelMixin):
    """
    Base calling for calling against simulated histograms
    """
    simulations = models.ForeignKey(SimulationsByCycles)

    class Meta:
        abstract = True

    @property
    def distance_metric(self):
        return pop_dist_corr_numpy

    @property
    def sim_hists_space(self):
        raise NotImplemented

    def find_best_in_space(self, hist):
        return get_closest(hist, self.sim_hists_space, self.distance_metric)

    def call_ms_hist(self, dbhist, microsatellite):
        return call_microsatellite_histogram(self, dbhist, microsatellite)


class FullMonoSimCorScheme(BaseSimCallingScheme, AllelesCyclesRangeMixin):
    """
    Calling schema for calling against simulated histograms
    """

    @property
    def sim_hists_space(self):
        yield from mono_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles)


class BaseBiAllelicMixin(object):

    @property
    def allele_number(self):
        return 2


class FullBiSimCorScheme(BaseSimCallingScheme, BaseBiAllelicMixin, FullRangeBiMixin):
    """
    Calling schema for calling against combinations of two simulated histograms
    """

    @property
    def sim_hists_space(self):
        yield from bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles)


class ProportionalSimCorSchemeMixin(object):

    @property
    def sim_hists_space(self):
        yield from proportional_bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles
        )


class ProportionalSimCorScheme(BestCorrelationProportionalCalledAlleleMixin, BaseBiAllelicMixin,
                               ProportionStepModelMixin, ProportionalAllelesCyclesRangeMixin, ProportionalSimCorSchemeMixin,
                               BaseSimCallingScheme):
    """
    Calling schema for calling against multi-allelic simulated histograms at differential proportions
    """

    pass  # the code is implemented in ProportionalSimCorSchemeMixin


class BoundProportionalSimCorScheme(BestCorrelationProportionalCalledAlleleMixin, ProportionStepModelMixin,
                                    ProportionsBoundsModelMixin, BaseBiAllelicMixin,
                                    BoundProportionalAllelesCyclesRangeMixin, ProportionalSimCorSchemeMixin,
                                    BaseSimCallingScheme):
    """
    Calling schema for calling against multi-allelic simulated histograms at differential proportions
    """

    pass  # the code is implemented in ProportionalSimCorSchemeMixin


class HighestPeaksBiSimCorSchemeModel(BestCorrelationProportionalHighestPeakCalledAlleleMIxin, ProportionStepModelMixin,
                                      ProportionsBoundsModelMixin, HighestPeaksRangeModelMixin, BaseBiAllelicMixin,
                                      BaseSimCallingScheme, FilterByHistMixin,
                                      BoundProportionalAllelesCyclesRangeMixin):

    @property
    def sim_hists_space(self):
        yield from proportional_bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles
        )
