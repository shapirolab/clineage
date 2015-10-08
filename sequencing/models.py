import os

from django.db import models
from django.conf import settings
from django.dispatch import receiver
from django.db.models.signals import post_save
from django.contrib.auth.models import User

from lib_prep.models import Library
from linapp.models import Protocol, RawData

### -------------------------------------------------------------------------------------
### Sequencing
### -------------------------------------------------------------------------------------
class MachineType(models.Model):
    company = models.CharField(max_length=50)
    model = models.CharField(max_length=50)
    def __unicode__(self):
        return self.company + '_' + self.model
### -------------------------------------------------------------------------------------
class Machine(models.Model):
    machineid = models.CharField(max_length=50)
    name = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(MachineType)
    ip = models.IPAddressField(null=True, blank=True)

    def __unicode__(self):
        return self.type.__unicode__() + '_' + self.machineid
### -------------------------------------------------------------------------------------
class NGSRun(models.Model):
    library = models.ForeignKey(Library)
    # data = models.ForeignKey(RawData, related_name='sequencing_event', null=True, blank=True)
    name = models.CharField(max_length=100, unique=True)
    machine = models.ForeignKey(Machine)
    protocol = models.ForeignKey(Protocol)
    user = models.ForeignKey(User)
    date = models.DateField()

#    def experiment(self):
#        experiments = []
#        for sample in self.samples.all():
#            if sample.experiment() not in experiments:
#                experiments.append(sample.experiment())
#        if len(experiments) <> 1:
#            raise
#        return experiments[0]
    def directory(self):
        path = '%s/%s'%(settings.NGS_RUNS,self.name)
        if not os.path.exists(path):
            os.makedirs(path)
        return path
    def __unicode__(self):
        return self.name
@receiver(post_save, sender=Sequencing)
def on_save(sender, **kwargs):
    kwargs['instance'].directory()
### -------------------------------------------------------------------------------------

class SequencingData(models.Model): # This contains the actual data.
    pass