from sequencing.calling.models import CallingScheme
from sequencing.calling.models_common import CalledAlleles
from django.db import models
import pickle


class CyclesMixin(models.Model):
    min_cycles = models.PositiveSmallIntegerField()
    max_cycles = models.PositiveSmallIntegerField()

    class Meta:
        abstract = True


class SimultaionsByCycles(CyclesMixin):
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


class SimCorScheme(CyclesMixin):
    simulations = models.ForeignKey(SimultaionsByCycles)


class BiallelicSimCorScheme(SimCorScheme):
    minimal_seeds_distance = models.PositiveSmallIntegerField()


class BestCor(CalledAlleles):
    confidence = models.FloatField()
    cycle = models.PositiveSmallIntegerField()
