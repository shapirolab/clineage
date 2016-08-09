from model_utils.managers import InheritanceManager
from sequencing.analysis.models import Histogram, MicrosatelliteHistogramGenotype
from django.db import models


class CallingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()

    def __str__(self):
        return self.name


class CalledGenotype(models.Model):
    histogram = models.ForeignKey(Histogram)
    microsatellite_genotype = models.ForeignKey(MicrosatelliteHistogramGenotype)
    calling_scheme = models.ForeignKey(CallingScheme)

    objects = InheritanceManager()


class HighestPeak(CalledGenotype):
    confidence = models.FloatField()