from sequencing.calling.models_common import CalledAlleles
from django.db import models


class HighestPeak(CalledAlleles):
    confidence = models.FloatField()