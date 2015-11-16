import os

from django.db import models
from django.contrib.auth.models import User
from django.dispatch import receiver
from django.db.models.signals import post_save
from django.core.urlresolvers import reverse
from django.conf import settings

from mptt.models import MPTTModel, TreeForeignKey
from utils.SequenceManipulations import *
#from genomes.models import Target, Sequence, TargetVariantType
#from lib_prep.models import Sequencing
#from wet_storage.models import *
#from sampling.models import Individual, Extraction, ExtractionEvent, SamplingEvent, CellContent, Cell
#from misc.models import *

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


### -------------------------------------------------------------------------------------
### Experiment Management
### -------------------------------------------------------------------------------------
# class Experiment(models.Model):
#    users = models.ManyToManyField(User, related_name='experiments', through='ExperimentUser')
#    name = models.CharField(max_length=50)
#    created_date = models.DateField(auto_now_add=True) #Automatically set the field to now when the object is first created.
#    description = models.TextField()
#    is_public = models.BooleanField(default=False)
#
#    def __unicode__(self):
#        return self.name
#    def individuals(self):
#        return Individual.objects.filter(
#            pk__in=list(set([extraction_event.individual.pk for extraction_event in self.extraction_events()]))
#        )
#    def extraction_events(self):
#        return ExtractionEvent.objects.filter(
#            pk__in=list(set([extraction.extraction_event.pk for extraction in self.extractions()]))
#        )
#    def extractions(self):
#        return Extraction.objects.filter(
#            pk__in=list(set([samplingevent.extraction.pk for samplingevent in self.samplingevents()]))
#        )
#    def samplingevents(self):
#        return SamplingEvent.objects.filter(
#            pk__in=list(set([cell.sampling.pk for cell in self.cells.all()]))
#        )
#    def cellscontents(self):
#        return CellContent.objects.filter(cell__in=self.cells.all())
#    # def sequencings(self):
#    #     return Sequencing.objects.filter(samples__in=self.cellscontents())
#    # def rawdata(self):
#    #     return RawData.objects.filter(sequencing__in=self.sequencings())
#    # def crd(self):
#    #     return CorrectedRawData.objects.filter(rawdata__in=self.rawdata())
#    def targetsanalysis(self):
#        return TargetAnalysis.objects.filter(sequencingdata__in=self.crd())
#    def targetsdistributions(self):  # Warning: returns list instead of queryset. TODO: reconsider that.
#        return [ta.distribution for ta in self.targetsanalysis()]
#    # def truesequences(self):
#    #     return TrueSequence.objects.filter(distribution__in=self.targetsdistributions())
#    def targetsvariants(self):
#        return TargetVariant.objects.filter(distribution__in=self.targetsdistributions())
#    def gensigs(self):
#        return GenSig.objects.filter(variants__in=self.targetsvariants())
#    def dms(self):
#        return DM.objects.filter(cell1__in=self.gensigs()).filter(cell2__in=self.gensigs())
#    def calculation(self):
#        return [crd.creation for crd in self.crd()].extend(
#            [ta.creation for ta in self.targetsanalysis()])
#        [crd.creation for crd in self.crd()]

### -------------------------------------------------------------------------------------
# class ExperimentLog(models.Model):
#    user = models.ForeignKey(User)
#    experiment = models.ForeignKey(Experiment, related_name='comments')
#    date = models.DateTimeField(auto_now=True)
#    comment = models.TextField()
### -------------------------------------------------------------------------------------
class FileContext(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name

### -------------------------------------------------------------------------------------
#def expriment_file_path(self, filename):
    #return settings.MEDIA_ROOT+'/CLineageFiles/'+self.experiment.name+'/'+filename
### -------------------------------------------------------------------------------------
# class ExperimentFile(models.Model):
#    def path(self, filename):
#        return settings.MEDIA_ROOT+'/CLineageFiles/'+self.experiment.name+'/'+filename
#    title = models.CharField(max_length=50)
#    experiment = models.ForeignKey(Experiment)
#    file_name = models.CharField(max_length=50)
#    file = models.FileField(upload_to=path)
#    context = models.ForeignKey(FileContext)
#    upload_date = models.DateField()
#    user = models.ForeignKey(User)
#    description = models.TextField()
#    def __unicode__(self):
#        return self.file_name
# ### -------------------------------------------------------------------------------------
# class ExperimentUser(models.Model):
#    user = models.ForeignKey(User)
#    experiment = models.ForeignKey(Experiment)
#    role = models.ForeignKey(LineageRole)
### -------------------------------------------------------------------------------------

class ProtocolType(models.Model):
    name = models.CharField(max_length=100)
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------

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
### -------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
#class Kit(models.Model):
    #TODO: complete.
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Algorithms description
### -------------------------------------------------------------------------------------
#class AlgorithmType(models.Model):#e.g. raw data to consensus data
    #name = models.CharField(max_length=50)
    #input = models.CharField(max_length=50)
    #output = models.CharField(max_length=50)
    #def __unicode__(self):
        #return self.name
#### -------------------------------------------------------------------------------------
#class Algorithm(models.Model):
    #name = models.CharField(max_length=50)
    #type = models.ForeignKey(AlgorithmType)
    #version = models.CharField(max_length=50)
    #developers = models.ManyToManyField(User)
    #def __unicode__(self):
        #return self.name

    #def get_absolute_url(self):
        #return reverse('algorithm-details', kwargs={'pk': self.pk})
#### -------------------------------------------------------------------------------------
#class AlgorithmParameter(models.Model):
    #algorithm = models.ManyToManyField(Algorithm, related_name='parameters')
    #name = models.CharField(max_length=50)
    #type = models.IntegerField() # todo: choices
    #def __unicode__(self):
        #return self.name
#### -------------------------------------------------------------------------------------
#class AlgorithmRun(models.Model):  # invocation
    #algorithm = models.ForeignKey(Algorithm, related_name='runs')
    #runname = models.CharField(max_length=50, null=True)
    #parameters = models.ManyToManyField(AlgorithmParameter, through='AlgorithmRunParameters')
    #user = models.ForeignKey(User)
    #timestamp = models.DateTimeField()
    #status = models.IntegerField() # todo: choices
    ## ExtraFiles = models.ManyToManyField(ExperimentFile)
    #def __unicode__(self):
        #return self.runname
    #def experiment(self):
##        if self.algorithm.type == AlgorithmType.objects.get() #TODO: Consider constraining the following 'if's according to AlgorithmType
        #if self.crd.all():
            #return self.crd.all()[0].experiment()
        #if self.targetsdistributions.all():
            #return self.targetsdistributions.all()[0].experiment()
        ## if self.truesequences.all():
        ##     return self.truesequences.all()[0].experiment()
        #if self.targetsvariants.all():
            #return self.targetsvariants.all()[0].experiment()
        #if self.gensigs.all():
            #return self.gensigs.all()[0].experiment()
        #if self.dms.all():
            #return self.dms.all()[0].experiment()
        #print 'ERROR: no experiment associated with SequenceDistribution!'
        #raise #TODO: implement proper exception here


#### -------------------------------------------------------------------------------------
#class AlgorithmRunParameters(models.Model):
    #run = models.ForeignKey(AlgorithmRun)
    #parameter = models.ForeignKey(AlgorithmParameter)
    #value = models.CharField(max_length=50)
### -------------------------------------------------------------------------------------
####### TODO: Log #############







#### -------------------------------------------------------------------------------------
#### Sequencing Data classes
#### -------------------------------------------------------------------------------------
#class RawData(models.Model):
    #sequencing = models.ForeignKey(Sequencing)
    #file = models.FilePathField(blank=True, null=True)
    #user = models.ForeignKey(User)
##    def __init__(self, *args, **kwargs):
##        super(RawData, self).__init__(*args, **kwargs)
##        if self.id: #non-empty instance
##            self.file.path = self.directory()
    #def __unicode__(self):
        #return '%s_raw'%unicode(self.sequencing)

    #def directory(self):
        #path = '%s/fastq'%(self.sequencing.directory())
        #if not os.path.exists(path):
            #os.makedirs(path)
        #return path
#### -------------------------------------------------------------------------------------
#class CorrectedRawData(models.Model):
    #rawdata = models.ForeignKey(RawData)
    #file = models.FilePathField()
    #creation = models.ForeignKey(AlgorithmRun, related_name = 'crd')

    #def __unicode__(self):
        #return '%s_corrected'%unicode(self.rawdata)

    #def experiment(self):
        #return self.rawdata.experiment()

    #def directory(self):
        #path = '%s/corrected'%(self.rawdata.sequencing.directory())
        #if not os.path.exists(path):
            #os.makedirs(path)
        #return path

#@receiver(post_save, sender=CorrectedRawData)
#def on_save(sender, **kwargs):
    #kwargs['instance'].directory()
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Targets Hierarchy
### -------------------------------------------------------------------------------------
#class FailedTargetValue(models.Model):
    #comment = models.TextField()
#### -------------------------------------------------------------------------------------
#class SequenceDistribution(models.Model):
    #target = models.ManyToManyField(Target, through='TargetAnalysis')
    #value = models.ForeignKey(Sequence)
    #count = models.IntegerField()
    #failed = models.ForeignKey(FailedTargetValue, null=True)

    #def __unicode__(self):
        #return self.value+'-'+str(self.count)

    #def experiment(self):
        #experiments = []
        #for analysis in self.targetanalysis_set.all():
            #experiments.append(analysis.sequencingdata.experiment())
        #experiment = list(set(experiments))
        #if len(experiment) > 1:
            #print 'ERROR: multiple experiments associated with SequenceDistribution!'
            #raise #TODO: implement proper exception here
        #if experiment:
            #return experiment[0]
        #print 'ERROR: no experiment associated with SequenceDistribution!'
        #raise #TODO: implement proper exception here
### -------------------------------------------------------------------------------------
#class TargetAnalysis(models.Model):
    #target = models.ForeignKey(Target)
    #distribution = models.ForeignKey(SequenceDistribution)
    #creation = models.ForeignKey(AlgorithmRun, related_name='targetsdistributions')
    #sequencingdata = models.ForeignKey(CorrectedRawData)

    #def __unicode__(self):
        #return self.target.__unicode__()+'_'+self.distribution.__unicode__()

    #def experiment(self):
        #return self.sequencingdata.experiment()
### -------------------------------------------------------------------------------------
# class TrueSequence(models.Model):
#     distribution = models.ManyToManyField(SequenceDistribution)
#     creation = models.ForeignKey(AlgorithmRun, related_name='truesequences')
#     value = models.ForeignKey(Sequence)
#
#     def __unicode__(self):
#         return self.value
#
#     def experiment(self):
#         experiments = []
#         for seqdist in self.distribution:
#             experiments.append(seqdist.experiment())
#         experiment = list(set(experiments))
#         if len(experiment) > 1:
#             print 'ERROR: multiple experiments associated with TrueSequence!'
#             raise #TODO: implement proper exception here
#         if experiment:
#             return experiment[0]
#         print 'ERROR: no experiment associated with TrueSequence!'
#         raise #TODO: implement proper exception here
### -------------------------------------------------------------------------------------
#class TargetVariant(models.Model):
    #distribution = models.ManyToManyField(SequenceDistribution)
    #creation = models.ForeignKey(AlgorithmRun, related_name='targetsvariants')
    #type = models.ForeignKey(TargetVariantType)
    #value = models.PositiveIntegerField()
    ##TODO:Decide if the same target as parent or different one.

    #def __unicode__(self):
        #return self.type.__unicode__()+'_'+str(self.value)

    #def experiment(self):
        #return self.ts.experiment()
### -------------------------------------------------------------------------------------
#class GenSig(models.Model):
    #variants = models.ManyToManyField(TargetVariant)
    #creation = models.ForeignKey(AlgorithmRun, related_name='gensigs')
    #value = models.TextField() #vector

    #def __unicode__(self):
        #return self.value

    #def experiment(self):
        #experiments = []
        #for tvariant in self.variants:
            #experiments.append(tvariant.experiment())
        #experiment = list(set(experiments))
        #if len(experiment) > 1:
            #print 'ERROR: multiple experiments associated with GenSig!'
            #raise #TODO: implement proper exception here
        #if experiment:
            #return experiment[0]
        #print 'ERROR: no experiment associated with GenSig!'
        #raise #TODO: implement proper exception here
#### -------------------------------------------------------------------------------------
#class DM(models.Model):
    #cell1 = models.ForeignKey(GenSig, related_name='+')
    #cell2 = models.ForeignKey(GenSig, related_name='+')
    #distance = models.IntegerField()
    #creation = models.ForeignKey(AlgorithmRun, related_name='dms')

    #def __unicode__(self):
        #return self.cell1.__unicode__()+'_'+self.cell2.__unicode__()+'='+str(self.distance)

    #def experiment(self):
        #if self.cell1.experiment() == self.cell2.experiment():
            #return self.cell1.experiment()
        #else:
            #print 'ERROR: multiple experiments associated with DM!'
            #raise #TODO: implement proper exception here
#### -------------------------------------------------------------------------------------
#class CellTreeNode(MPTTModel):
    #cell = models.ForeignKey(Cell)
    #parent = TreeForeignKey('self', null=True, blank = True, related_name='children')
    #distance = models.IntegerField()
    #def __unicode__(self):
        #return self.cell.name
##    class MPTTMeta:
##        order_insertion_by=['name']
### -------------------------------------------------------------------------------------



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
