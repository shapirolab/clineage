from django.db import models
from sequencing.calling.simcor.simulation_spaces import get_far_apart_highest_peaks


class TrimmedSeedsSearchRangeMixin(object):
    def trim_seeds(self, hist):
        get_far_apart_highest_peaks(hist=hist)


class TrimmedSeedsBiallelicSearchRangeMixin(models.Model, TrimmedSeedsSearchRangeMixin):
    minimal_seeds_distance = models.PositiveSmallIntegerField()

    class Meta:
        abstract = True

    def trim_seeds(self, hist):
        get_far_apart_highest_peaks(hist=hist, k=2, d=self.minimal_seeds_distance)
