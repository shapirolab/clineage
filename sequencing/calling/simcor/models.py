import pickle
import itertools
from django.db import models
from sequencing.calling.models import CallingScheme
from sequencing.calling.models_common import CalledAlleles
from sequencing.calling.simcor.simulation_spaces import mono_sim_hists_space_generator, bi_sim_hists_space_generator,\
    proportional_bi_sim_hists_space_generator
from sequencing.calling.simcor.seeds_range_trimming import TrimmedSeedsBiallelicSearchRangeMixin

class CyclesMixin(models.Model):
    min_cycles = models.PositiveSmallIntegerField()
    max_cycles = models.PositiveSmallIntegerField()

    @property
    def cycle_bounds(self):
        return self.min_cycles, self.max_cycles

    class Meta:
        abstract = True


class SimulationsByCycles(CyclesMixin):
    """
    A simulations pickle file
     A dictionary for a specific MS type, indexed like:
      d[ms_len][cycle]
    """
    name = models.CharField(max_length=50)
    sim_hists = models.FilePathField(max_length=200)
    description = models.TextField()

    def __str__(self):
        return self.name

    def get_simulations_dict(self):
        with open(self.sim_hists, 'rb') as f:
            data = pickle.load(f)
            return data


class SimCorScheme(CallingScheme, CyclesMixin):
    """
    Base calling schema for calling against simulated histograms
    """
    simulations = models.ForeignKey(SimulationsByCycles)

    @property
    def seeds_and_cycles(self):
        yield from itertools.product(range(*self.ms_len_bounds), range(*self.cycle_bounds))

    @property
    def sim_hists_space(self):
        yield from mono_sim_hists_space_generator(
            self.simulations,
            self.seeds_and_cycles)


class NaiveBiallelicSimCorScheme(SimCorScheme, TrimmedSeedsBiallelicSearchRangeMixin):
    """
    Biallelic calling schema
     Taking into account a trimmed simulation space based on highest peaks
    """

    @property
    def sim_hists_space(self):
        yield from bi_sim_hists_space_generator(
            self.simulations,
            self.seeds_and_cycles)


# class NaiveBiallelicSimCorScheme(SimCorScheme, TrimmedSeedsBiallelicSearchRangeMixin):
#     """
#     Biallelic calling schema
#      Taking into account a trimmed simulation space based on highest peaks
#     """
#
#     @property
#     def sim_hists_space(self):
#         yield from bi_sim_hists_space_generator(
#             self.simulations,
#             self.seeds_and_cycles)


class BestCorrelationCalledAlleles(CalledAlleles):
    """
    Calling result with confidence (1-correlation) and simulated cycle.
    """
    confidence = models.FloatField()
    cycle = models.PositiveSmallIntegerField()
