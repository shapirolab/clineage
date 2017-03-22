import itertools
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


class ProportionStepModelMixin(models.Model):
    proportion_step = models.DecimalField(max_digits=3, decimal_places=2)

    class Meta:
        abstract = True


class ProportionsBoundsModelMixin(models.Model, ProportionsBoundsMixin):
    lower_prop_bound = models.DecimalField(max_digits=3, decimal_places=2)
    upper_prop_bound = models.DecimalField(max_digits=3, decimal_places=2)

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
        return cls.get_for_alleles([allele])

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


class Proportions(models.Model):
    p1 = models.DecimalField(max_digits=3, decimal_places=2, null=True)
    p2 = models.DecimalField(max_digits=3, decimal_places=2, null=True)
    p3 = models.DecimalField(max_digits=3, decimal_places=2, null=True)
    p4 = models.DecimalField(max_digits=3, decimal_places=2, null=True)
    _num_of_allele_fields = 4

    @classmethod
    def proportion_field_names(cls):
        for i in range(1, cls._num_of_allele_fields + 1):
            yield 'p{}'.format(i)

    @classmethod
    def get_for_proportions(cls, proportions):
        assert len(list(proportions)) <= cls._num_of_allele_fields
        names = list(cls.proportion_field_names())
        assert len(names) >= len(proportions)
        ordered_proportions = dict(itertools.zip_longest(
            names,
            proportions,
            fillvalue=None,
        ))
        obj, c = cls.objects.get_or_create(**ordered_proportions)
        return obj

    @property
    def proportion_fields(self):
        return [
            self.p1,
            self.p2,
            self.p3,
            self.p4,
        ]

    @property
    def proportions(self):
        return [p for p in self.proportion_fields if p is not None]


class ProportionalMicrosatelliteAlleleSet(MicrosatelliteAlleleSet):
    proportions = models.ForeignKey(Proportions)

    @classmethod
    def get_for_proportional_alllels(cls, mas, proportions_dict):
        pmas = cls(microsatellitealleleset_ptr_id=mas.pk)
        pmas.__dict__.update(mas.__dict__)
        ordered_proportions = [proportions_dict[a] for a in mas.alleles]
        pmas.proportions = Proportions.get_for_proportions(ordered_proportions)
        pmas.save()
        return pmas

    @property
    def alleles(self):
        return frozenset(zip([a for a in self.allele_fields if a is not None], self.proportions.proportions))


class BestCorrelationProportionalCalledAlleles(BestCorrelationCalledAlleles):
    """
    Calling result with confidence (1-correlation) and simulated cycle.
    """
    pass


class BestCorrelationProportionalHighestPeakCalledAlleles(BestCorrelationCalledAlleles):
    """
    Calling result with confidence (1-correlation) and simulated cycle.
    """
    pass
