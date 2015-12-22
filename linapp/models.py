import os

from django.db import models
from django.contrib.auth.models import User
from django.dispatch import receiver
from django.db.models.signals import post_save
from django.core.urlresolvers import reverse
from django.conf import settings

from mptt.models import MPTTModel, TreeForeignKey
from utils.SequenceManipulations import *
#from sampling.models import Individual, Cell

# Create your models here.

### -------------------------------------------------------------------------------------
### Users/Roles Management
### -------------------------------------------------------------------------------------
class LineageRole(models.Model):
    name = models.CharField(max_length=50)
    read = models.BooleanField(default=True)
    write = models.BooleanField(default=False)
    delete = models.BooleanField(default=False)
    def __unicode__(self):
        return self.name
    @staticmethod
    def emptyRole():
        return LineageRole(name = 'Empty', read = False, write = False, delete = False)
    @staticmethod
    def emptySuperUserRole():
        return LineageRole(name = 'Empty', read = True, write = True, delete = True)
### -------------------------------------------------------------------------------------
class UserProfile(models.Model):
    user = models.OneToOneField(User)
    institute = models.CharField(max_length=50)
    comment = models.TextField()
    def create_user_profile(sender, instance, created, **kwargs):
        if created:
            UserProfile.objects.create(user=instance)

    post_save.connect(create_user_profile, sender=User)
### -------------------------------------------------------------------------------------



class ProtocolType(models.Model):
    name = models.CharField(max_length=100)
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Protocol(models.Model):
    initials = models.CharField(max_length=10)
    name = models.CharField(max_length=100)
    abstract = models.TextField()
    fulldescription = models.TextField()
    kit = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(ProtocolType)
    file = models.FilePathField(null=True, blank=True)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------

# NOTE: these are required by linapp.migrations.0001_initial
def expriment_file_path(self, filename):
    return settings.MEDIA_ROOT+'/CLineageFiles/'+self.experiment.name+'/'+filename

def sampling_event_path(self, filename):
    return settings.MEDIA_ROOT+'/CLineageFiles/sampling_events/'+self.name+'/'+filename

###----------------------------------------------------------------------------------------
### User Report Data and Tables
###----------------------------------------------------------------------------------------
class UserReport(models.Model):
    cells = models.ManyToManyField('sampling.Cell')
    partner = models.ForeignKey(User, null=True)
    individual = models.ManyToManyField('sampling.Individual', null=True)
    creation_date = models.DateField(auto_now_add=True) #Automatically set the field to now when the object is first created.

    def __unicode__(self):
        return self.partner.username if self.partner else 'detached'

    @staticmethod
    def get_create_new(cells, partner=None, individual=None):
        user_report = None
        for ur in UserReport.objects.all():
            if set(ur.cells.all()) == set(cells):
                user_report = ur
        if user_report:
            print user_report
            if not user_report.partner and partner:
                user_report.partner = partner
                user_report.save()
            if not user_report.individual and individual:
                user_report.individual = individual
                user_report.save()
            return user_report
        else:
            if not cells:
                raise
            ur = UserReport.objects.create(partner=partner)
            ur.individual = individual
            ur.cells = cells
            return ur

class Sequence(models.Model):
    def create(_length=0, _sequence='', _hash=hashlib.md5('').hexdigest()):
        if _length > 0 and _sequence != '' and _hash != ''\
           and _length == len(_sequence) and _hash == hashlib.md5(_sequence).hexdigest():
            return Sequence(length=_length, sequence=_sequence, hash=_hash)
        return Sequence(length=len(_sequence), sequence=_sequence, hash=hashlib.md5(_sequence).hexdigest())
    create = staticmethod(create)

    length = models.IntegerField()
    sequence = models.TextField()
    hash = models.CharField(max_length=32, unique=True) #md5(sequence) for enabling uniqueness and fast comparison.

