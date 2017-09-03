from model_utils.managers import InheritanceManager
from sequencing.analysis.models import Histogram
from targeted_enrichment.planning.models import Microsatellite
from django.db import models
import itertools
from sequencing.calling.range import MultiAlleleMixin


class CallingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()

    objects = InheritanceManager()

    def __str__(self):
        return self.name

    @property
    def called_allele_class(self):
        raise NotImplemented

    @property
    def called_allele_class(self):
        raise NotImplemented

    @property
    def call_ms_hist(self, dbhist, microsatellite):
        raise NotImplemented


class MicrosatelliteAlleleSet(models.Model, MultiAlleleMixin):
    allele1 = models.PositiveSmallIntegerField(default=0)  # use 0 instead of null to allow unique enforcement by MySQL
    allele2 = models.PositiveSmallIntegerField(default=0)  # use 0 instead of null to allow unique enforcement by MySQL
    allele3 = models.PositiveSmallIntegerField(default=0)  # use 0 instead of null to allow unique enforcement by MySQL
    allele4 = models.PositiveSmallIntegerField(default=0)  # use 0 instead of null to allow unique enforcement by MySQL
    _num_of_allele_fields = 4

    objects = InheritanceManager()

    @staticmethod
    def sort_alleles(alleles):
        return sorted(list(alleles))

    @property
    def allele_number(self):
        return self._num_of_allele_fields

    @classmethod
    def allele_field_names(cls):
        for i in range(1, cls._num_of_allele_fields+1):
            yield 'allele{}'.format(i)

    @classmethod
    def get_for_alleles(cls, alleles):
        assert len(alleles) <= cls._num_of_allele_fields
        l = cls.sort_alleles(alleles)
        names = list(cls.allele_field_names())
        assert len(names) >= len(l)
        ordered_genotypes = dict(itertools.zip_longest(
            names,
            l,
            fillvalue=0,
        ))
        obj, c = cls.objects.get_or_create(**ordered_genotypes)
        return obj

    @property
    def allele_fields(self):
        return [
            self.allele1,
            self.allele2,
            self.allele3,
            self.allele4,
        ]

    @property
    def alleles(self):
        return {
            allele for allele
            in self.allele_fields
            if allele != 0
        }

    class Meta:
        unique_together = (
            (
                "allele1",
                "allele2",
                "allele3",
                "allele4",
            ),
        )


class CalledAlleles(models.Model):
    histogram = models.ForeignKey(Histogram)
    microsatellite = models.ForeignKey(Microsatellite)
    genotypes = models.ForeignKey(MicrosatelliteAlleleSet)
    calling_scheme = models.ForeignKey(CallingScheme)

    objects = InheritanceManager()

    class Meta:
        unique_together = (
            (
                "histogram",
                "microsatellite",
                "calling_scheme",
            ),
        )


