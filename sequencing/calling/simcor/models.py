from sequencing.calling.models import CallingScheme
from sequencing.calling.models_common import CalledAlleles
from django.db import models
import pickle


class SimultaionsByCycles(models.Model):
    """
    A simulations pickle file
     A dictionary for a specific MS type, indexed like:
      d[ms_len][cycle]
    """
    name = models.CharField(max_length=50)
    min_cycles = models.PositiveIntegerField()
    max_cycles = models.PositiveIntegerField()
    sim_hists = models.FilePathField(max_length=200)
    description = models.TextField()

    def __str__(self):
        return self.name

    def get_simulations_dict(self):
        with open(self.sim_hists, 'rb') as f:
            data = pickle.load(f)
            return data


class BestCor(CalledAlleles):
    confidence = models.FloatField()
    cycle = models.PositiveSmallIntegerField()
    simulations = models.ForeignKey(SimultaionsByCycles)
