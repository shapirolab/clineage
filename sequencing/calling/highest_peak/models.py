from sequencing.calling.models_common import CalledGenotype
from django.db import models


class HighestPeak(CalledGenotype):
    confidence = models.FloatField()