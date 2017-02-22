from model_utils.managers import InheritanceManager
from sequencing.analysis.models import Histogram
from targeted_enrichment.planning.models import Microsatellite
from django.db import models
import itertools


class CallingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()

    def __str__(self):
        return self.name


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


class MicrosatelliteAlleleSet(models.Model):
    allele1 = models.PositiveSmallIntegerField()
    allele2 = models.PositiveSmallIntegerField()
    allele3 = models.PositiveSmallIntegerField()
    allele4 = models.PositiveSmallIntegerField()
    _num_of_allele_fields = 4

    objects = InheritanceManager()

    @classmethod
    def allele_field_names(cls):
        for i in range(1, cls._num_of_allele_fields+1):
            yield 'microsatellite_genotype{}'.format(i)

    @classmethod
    def get_for_repeats(cls, alleles):
        assert len(alleles) <= cls._num_of_allele_fields
        l = sorted(list(alleles))
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


class SingleMicrosatelliteAlleleSet(MicrosatelliteAlleleSet):
    _num_of_allele_fields = 1

    @classmethod
    def get_for_repeat(cls, allele):  # Optional, TODO: consider removal
        assert len(allele) == SingleMicrosatelliteAlleleSet._num_of_allele_fields
        return cls.get_for_repeats([allele])

    @property
    def allele_fields(self):
        return [
            self.allele1,
        ]


class ProportionalMSAlleleSet(MicrosatelliteAlleleSet):
    p1 = models.DecimalField(max_digits=3, decimal_places=2)
    p2 = models.DecimalField(max_digits=3, decimal_places=2)
    p3 = models.DecimalField(max_digits=3, decimal_places=2)
    p4 = models.DecimalField(max_digits=3, decimal_places=2)
    _num_of_allele_fields = 4

    @property
    def proportion_fields(self):
        return [
            self.p1,
            self.p2,
            self.p3,
            self.p4,
        ]

    @property
    def alleles(self):
        return {
            zip(
                [allele for allele
                in self.allele_fields
                if allele != 0],
                [p for p
                 in self.proportion_fields],
            )
        }


class ProportionalBiAllelicMSAlleleSet(ProportionalMSAlleleSet):
    _num_of_allele_fields = 2

    @property
    def proportion_fields(self):
        return [
            self.p1,
            self.p2,
        ]