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
    HighestPeaksRangeMixin
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

    def filter_by_hist(self, hist):
        raise NotImplemented

    def find_best_in_space(self, hist):
        return get_closest(hist, self.filter_by_hist(hist, self.sim_hists_space), self.distance_metric)


class FilterByHistMixin(DynamicFilteredHistSpaceMixin):

    @property
    def filter_by_hist(self, hist):
        """cuts the simulations based on the hist"""
        alleles_by_hist = self.alleles_by_hist(hist)
        yield from filterfalse(lambda self: (allele for allele in self.sim_hists_space.allele_frozenset in alleles_by_hist),
                               self.sim_hists_space)


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


class ProportionalSimCorScheme(BestCorrelationProportionalCalledAlleleMixin, BaseBiAllelicMixin, BaseSimCallingScheme,
                               ProportionStepModelMixin, ProportionalAllelesCyclesRangeMixin):
    """
    Calling schema for calling against multi-allelic simulated histograms at differential proportions
    """

    @property
    def sim_hists_space(self):
        yield from proportional_bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles
        )


class BoundProportionalSimCorScheme(BestCorrelationProportionalCalledAlleleMixin, ProportionStepModelMixin,
                                    ProportionsBoundsModelMixin, BaseBiAllelicMixin, BaseSimCallingScheme,
                                    BoundProportionalAllelesCyclesRangeMixin):
    """
    Calling schema for calling against multi-allelic simulated histograms at differential proportions
    """

    @property
    def sim_hists_space(self):
        yield from proportional_bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles
        )


class HighestPeaksBiSimCorScheme(
                                 BestCorrelationProportionalHighestPeakCalledAlleleMIxin,
                                 ProportionStepModelMixin,
                                 ProportionsBoundsModelMixin,
                                 FilterByHistMixin,
                                 BaseBiAllelicMixin,
                                 BaseSimCallingScheme,
                                 BoundProportionalAllelesCyclesRangeMixin,
                                 HighestPeaksRangeMixin
                                ):

    @property
    def sim_hists_space(self):
        yield from proportional_bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles
        )




    # @property
    # def sim_hists_space(self):
    #     filterfalse(lambda x: x in self.alleles_by_hist, proportional_bi_sim_hists_space_generator(
    #         self.simulations.get_simulations_dict(),
    #         self.alleles_and_cycles
    #     ))



# class NaiveBiallelicSimCorScheme(BaseSimCorMixin, TrimmedSeedsBiallelicSearchRangeMixin):
#     """
#     Biallelic calling schema
#      Taking into account a trimmed simulation space based on highest peaks
#     """
#
#     @property
#     def sim_hists_space(self):
#         yield from bi_sim_hists_space_generator(
#             self.simulations.get_simulations_dict(),
#             self.seeds_and_cycles)
