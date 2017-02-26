import pickle
from django.db import models
from sequencing.calling.range import MSLengthBoundsMixin, ProportionsBoundsMixin
from sequencing.calling.models_common import MicrosatelliteAlleleSet
from sequencing.calling.models_common import CalledAlleles


class CycleBoundsMixin(object):

    @property
    def cycle_bounds(self):
        raise NotImplemented


class CyclesRangeMixin(CycleBoundsMixin):

    @property
    def cycles(self):
        yield from range(*self.cycle_bounds)


class CyclesModelMixin(models.Model, CycleBoundsMixin):
    min_cycles = models.PositiveSmallIntegerField()
    max_cycles = models.PositiveSmallIntegerField()

    @property
    def cycle_bounds(self):
        return self.min_cycles, self.max_cycles

    class Meta:
        abstract = True


class SimulationsByCycles(CyclesModelMixin):
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


class MSLengthBoundsModelMixin(models.Model, MSLengthBoundsMixin):
    min_ms_len = models.PositiveSmallIntegerField()
    max_ms_len = models.PositiveSmallIntegerField()

    @property
    def ms_len_bounds(self):
        return self.min_ms_len, self.max_ms_len

    class Meta:
        abstract = True


class ProportionsBoundsModelMixin(models.Model, ProportionsBoundsMixin):
    lower_prop_bound = models.DecimalField(max_digits=2, decimal_places=2)
    upper_prop_bound = models.DecimalField(max_digits=2, decimal_places=2)

    @property
    def proportion_bounds(self):
        return self.lower_prop_bound, self.upper_prop_bound

    class Meta:
        abstract = True


class SingleMicrosatelliteAlleleSet(MicrosatelliteAlleleSet):
    _num_of_allele_fields = 1

    @classmethod
    def get_for_repeat(cls, allele):
        assert len(allele) == cls._num_of_allele_fields
        return cls.get_for_repeats([allele])

    @property
    def allele_fields(self):
        return [
            self.allele1,
        ]


# class ProportionalBiAllelicMSAlleleSet(ProportionalMSAlleleSet):
#     _num_of_allele_fields = 2
#
#     @property
#     def proportion_fields(self):
#         return [
#             self.p1,
#             self.p2,
#         ]


class BestCorrelationCalledAlleles(CalledAlleles):
    """
    Calling result with confidence (1-correlation) and simulated cycle.
    """
    confidence = models.FloatField()
    cycle = models.PositiveSmallIntegerField()


class BestCorrelationProportionalCalledAlleles(BestCorrelationCalledAlleles):
    """
    Calling result with confidence (1-correlation) and simulated cycle.
    """
    pass
