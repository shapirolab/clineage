from django.db import models
from sequencing.calling.models import CallingScheme
from sequencing.calling.simcor.simulation_spaces import mono_sim_hists_space_generator, bi_sim_hists_space_generator,\
    proportional_bi_sim_hists_space_generator
from sequencing.calling.hist_dist import pop_dist_corr_numpy
from sequencing.calling.simcor.calling import call_microsatellite_histogram
from sequencing.calling.simcor.models_common import CyclesModelMixin, SimulationsByCycles, MSLengthBoundsModelMixin, \
    ProportionsBoundsModelMixin
from sequencing.calling.simcor.range import AllelesCyclesRangeMixin, FullRangeBiMixin, \
    ProportionalAllelesCyclesRangeMixin, BoundProportionalAllelesCyclesRangeMixin


class BaseSimCallingScheme(CallingScheme, MSLengthBoundsModelMixin, CyclesModelMixin):
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

    def call_ms_hist(self, dbhist, microsatellite):
        raise NotImplemented


class FullMonoSimCorScheme(BaseSimCallingScheme, AllelesCyclesRangeMixin):
    """
    Calling schema for calling against simulated histograms
    """

    @property
    def sim_hists_space(self):
        yield from mono_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles)

    def call_ms_hist(self, dbhist, microsatellite):
        return call_microsatellite_histogram(self, dbhist, microsatellite)


class FullBiSimCorScheme(BaseSimCallingScheme, FullRangeBiMixin):
    """
    Calling schema for calling against combinations of two simulated histograms
    """
    @property
    def allele_number(self):
        return 2

    @property
    def sim_hists_space(self):
        yield from bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles)

    def call_ms_hist(self, dbhist, microsatellite):
        return call_microsatellite_histogram(self, dbhist, microsatellite)


class ProportionalSimCorScheme(BaseSimCallingScheme,
                               ProportionalAllelesCyclesRangeMixin):
    """
    Calling schema for calling against multi-allelic simulated histograms at differential proportions
    """

    @property
    def sim_hists_space(self):
        yield from proportional_bi_sim_hists_space_generator(
            self.simulations.get_simulations_dict(),
            self.alleles_and_cycles
        )

    def call_ms_hist(self, dbhist, microsatellite):
        return call_microsatellite_histogram(dbhist, microsatellite, self)


class BoundProportionalSimCorScheme(BaseSimCallingScheme, ProportionsBoundsModelMixin,
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

    def call_ms_hist(self, dbhist, microsatellite):
        return call_microsatellite_histogram(dbhist, microsatellite, self)


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
