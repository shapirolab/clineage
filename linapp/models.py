import hashlib
import mmap
import os

from django.db import models, transaction
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic
from django.contrib.auth.models import User
from django.dispatch import receiver
from django.db.models import Count
from django.db.models.signals import *
from django.core.urlresolvers import reverse
from django.conf import settings
from mptt.models import MPTTModel, TreeForeignKey
from utils.wells import index2str, str2index
from utils.SequenceManipulations import *

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
class Experiment(models.Model):
    users = models.ManyToManyField(User, related_name='experiments', through='ExperimentUser')
    name = models.CharField(max_length=50)
    created_date = models.DateField(auto_now_add=True) #Automatically set the field to now when the object is first created.
    description = models.TextField()
    is_public = models.BooleanField(default=False)

    def __unicode__(self):
        return self.name
    def individuals(self):
        return Individual.objects.filter(
            pk__in=list(set([extraction_event.individual.pk for extraction_event in self.extraction_events()]))
        )
    def extraction_events(self):
        return ExtractionEvent.objects.filter(
            pk__in=list(set([extraction.extraction_event.pk for extraction in self.extractions()]))
        )
    def extractions(self):
        return Extraction.objects.filter(
            pk__in=list(set([samplingevent.extraction.pk for samplingevent in self.samplingevents()]))
        )
    def samplingevents(self):
        return SamplingEvent.objects.filter(
            pk__in=list(set([cell.sampling.pk for cell in self.cells.all()]))
        )
    def cellscontents(self):
        return CellContent.objects.filter(cell__in=self.cells.all())
    def sequencings(self):
        return Sequencing.objects.filter(samples__in=self.cellscontents())
    def rawdata(self):
        return RawData.objects.filter(sequencing__in=self.sequencings())
    def crd(self):
        return CorrectedRawData.objects.filter(rawdata__in=self.rawdata())
    def targetsanalysis(self):
        return TargetAnalysis.objects.filter(sequencingdata__in=self.crd())
    def targetsdistributions(self):  # Warning: returns list instead of queryset. TODO: reconsider that.
        return [ta.distribution for ta in self.targetsanalysis()]
    # def truesequences(self):
    #     return TrueSequence.objects.filter(distribution__in=self.targetsdistributions())
    def targetsvariants(self):
        return TargetVariant.objects.filter(distribution__in=self.targetsdistributions())
    def gensigs(self):
        return GenSig.objects.filter(variants__in=self.targetsvariants())
    def dms(self):
        return DM.objects.filter(cell1__in=self.gensigs()).filter(cell2__in=self.gensigs())
#    def calculation(self):
#        return [crd.creation for crd in self.crd()].extend(
#            [ta.creation for ta in self.targetsanalysis()])
#        [crd.creation for crd in self.crd()]

### -------------------------------------------------------------------------------------
class ExperimentLog(models.Model):
    user = models.ForeignKey(User)
    experiment = models.ForeignKey(Experiment, related_name='comments')
    date = models.DateTimeField(auto_now=True)
    comment = models.TextField()
### -------------------------------------------------------------------------------------
class FileContext(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class ExperimentFile(models.Model):
    def path(self, filename):
        return settings.MEDIA_ROOT+'/CLineageFiles/'+self.experiment.name+'/'+filename
    title = models.CharField(max_length=50)
    experiment = models.ForeignKey(Experiment)
    file_name = models.CharField(max_length=50)
    file = models.FileField(upload_to=path)
    context = models.ForeignKey(FileContext)
    upload_date = models.DateField()
    user = models.ForeignKey(User)
    description = models.TextField()
    def __unicode__(self):
        return self.file_name
### -------------------------------------------------------------------------------------
class ExperimentUser(models.Model):
    user = models.ForeignKey(User)
    experiment = models.ForeignKey(Experiment)
    role = models.ForeignKey(LineageRole)
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Types and descriptors
### -------------------------------------------------------------------------------------
class ProtocolType(models.Model):
    name = models.CharField(max_length=100)
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class SampleComposition(models.Model):#e.g. single cell or bulk
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class CellContentType(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class TargetType(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class TargetVariantType(models.Model): #TODO: This might be useless. Decide.
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
 ### -------------------------------------------------------------------------------------
class SampleStatus(models.Model):
    name = models.CharField(max_length=50)

    class Meta:
        verbose_name = 'Sample status'
        verbose_name_plural = 'Samples status'

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class TargetEnrichmentFailureType(models.Model):
    """
    1 = No product
    2 = Primer dimer is wider, equal or  close to the same band width of expected product
    3 = Smear or more than 3 products  (other than the primer dimer which can be purified). If less than 3, real product
        has to be wider than byproducts.
    4 = More than 1 product is in the range of correct size.
    5 = NGS failure - Primer pair did not work in the context of a successful NGS run (amplified and sequenced).
    """
    name = models.CharField(max_length=50)
    description = models.TextField(null=True, blank=True)
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Generic Biological Objects
### -------------------------------------------------------------------------------------

class Sequence(models.Model):
    def create(_length=0, _sequence='', _hash=hashlib.md5('').hexdigest()):
        if _length > 0 and _sequence <> '' and _hash <> ''\
           and _length == len(_sequence) and _hash == hashlib.md5(_sequence).hexdigest():
            return Sequence(length=_length, sequence=_sequence, hash=_hash)
        return Sequence(length=len(_sequence), sequence=_sequence, hash=hashlib.md5(_sequence).hexdigest())
    create = staticmethod(create)

    length = models.IntegerField()
    sequence = models.TextField()
    hash = models.CharField(max_length=32, unique=True) #md5(sequence) for enabling uniqueness and fast comparison.
### -------------------------------------------------------------------------------------
class Taxa(models.Model):
    name = models.CharField(max_length=50)
    taxonomy_id = models.IntegerField()
    rank = models.CharField(max_length=50)
    parent = models.IntegerField(null=True, blank=True)
    friendly_name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class GeneticBackground(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Organ(models.Model):
    name = models.CharField(max_length=50)
    #rank = models.IntegerField()
    #parent = models.IntegerField()

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Tissue(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Assembly(models.Model):
    taxa = models.ForeignKey(Taxa)
    name = models.CharField(max_length=50)
    friendly_name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.friendly_name

    class Meta:
        verbose_name_plural = 'Assemblies'

    def get_path(self):
        return os.path.join(self.taxa.friendly_name, self.friendly_name)

### -------------------------------------------------------------------------------------
class Chromosome(models.Model):
    name = models.CharField(max_length=50)
    assembly = models.ForeignKey(Assembly)
    sequence_length = models.IntegerField(null=True)
    cyclic = models.BooleanField()

    def __unicode__(self):
        return self.name

    def get_path(self, ext="txt"):
        return os.path.join(self.assembly.get_path(), 'chr{}.{}'.format(self.name, ext))

    def get_abs_path(self):
        return os.path.join(settings.CHROMOSOMES_PATH, self.get_path())

    def getdna(self, start, stop):
        if start <= stop and stop <= self.sequence_length:
            with open(self.get_abs_path(), 'r+b') as f:
                mm = mmap.mmap(f.fileno(), 0)
                return mm[start-1:stop].upper()
        if self.cyclic and start > stop:
            return self.getdna(start, self.sequence_length) + self.getdna(0, stop)
        if self.cyclic and stop > self.sequence_length:
            return self.getdna(start, stop-self.sequence_length)
        raise ValueError('indices out of bounds')

    def locate(self, start, stop, sequence, padding=10):
        l = padding if start > padding else start
        r = padding if stop + padding < self.sequence_length else self.sequence_length - stop
        s = self.getdna(start - l, stop + r)

        # check for exact match
        if s[l:-r] == sequence.upper():
            return start, stop

        index = s.find(sequence.upper())  # search for reference value with extra 10bp on each side
        if index < 0:
            raise ValueError('could not find %s within %s:%d-%d' % (sequence, self.name, start-10, stop+10))

        return start - l + index, start - l + index + len(sequence) - 1

### -------------------------------------------------------------------------------------
#class Kit(models.Model):
    #TODO: complete.
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
class Target(models.Model):#Target is a locus on a reference genome.
    name = models.CharField(max_length=50)
    type = models.ForeignKey(TargetType) #Microsatellite / SNP / etc...
    chromosome = models.ForeignKey(Chromosome)
    start_pos = models.IntegerField(db_index=True)
    end_pos = models.IntegerField(db_index=True)
    referencevalue = models.ForeignKey(Sequence)

    # unique_together = (("chromosome", "start_pos"),)

    def get_referencevalue(self):
        return self.chromosome.getdna(self.start_pos, self.end_pos)

    def __unicode__(self):
        return self.name

    def validate_reference(self):
        assert self.referencevalue.sequence == self.get_referencevalue()

    def update_primers(self):
        for te in TargetEnrichment.objects.filter(chromosome=self.chromosome)\
                                .filter(left__end_pos__lte=self.start_pos)\
                                .filter(right__start_pos__gte=self.end_pos):
            te.update_enriched_targets()
### -------------------------------------------------------------------------------------
class PrimerTail(models.Model):
    tail = models.CharField(max_length=50, null=True)
### -------------------------------------------------------------------------------------
class Primer(Target):
    PLUS = '+'
    MINUS = '-'
    STRANDS = (
        (PLUS, 'Plus'),
        (MINUS, 'Minus'),
    )
    strand = models.CharField(max_length=1, choices=STRANDS, null=True)
    sequence = models.ForeignKey(Sequence)
    tail = models.ForeignKey(PrimerTail, null=True)
    physical_locations = generic.GenericRelation('SampleLocation',
                                             content_type_field='content_type',
                                             object_id_field='object_id')
    def validate_reference(self):
        if self.strand == self.PLUS:
            assert self.referencevalue.sequence == self.get_referencevalue()
            assert self.sequence.sequence[-(self.end_pos-self.start_pos+1):] == self.get_referencevalue()
        if self.strand == self.MINUS:
            assert self.referencevalue.sequence == self.get_referencevalue()
            assert complement(self.sequence.sequence[-(self.end_pos-self.start_pos+1):])[::-1] == self.get_referencevalue()
### -------------------------------------------------------------------------------------
class Microsatellite(Target):
    repeat_unit_len = models.PositiveIntegerField() #length of repeat Nmer
    repeat_unit_type = models.CharField(max_length=50) #string of repeat Nmer
    repeat_number = models.DecimalField(max_digits=5, decimal_places=1, null=True)
    def __unicode__(self):
        return self.name

### -------------------------------------------------------------------------------------
class TargetEnrichmentType(models.Model):
    name = models.CharField(max_length=50)
    protocol = models.ForeignKey(Protocol, null=True)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class TargetEnrichment(models.Model):
    type = models.ForeignKey(TargetEnrichmentType)
    chromosome = models.ForeignKey(Chromosome)
    left = models.ForeignKey(Primer, related_name='left_primer')
    right = models.ForeignKey(Primer, related_name='right_primer')
    amplicon = models.CharField(max_length=500)
    passed_validation = models.NullBooleanField()
    validation_failure = models.ForeignKey(TargetEnrichmentFailureType, null=True)
    validation_date = models.DateField(null=True, blank=True)
    comment = models.CharField(max_length=50, blank=True, null=True)
    physical_locations = generic.GenericRelation('SampleLocation',
                                                 content_type_field='content_type',
                                                 object_id_field='object_id')
    targets = models.ManyToManyField(Target, related_name='primer_pair', null=True, blank=True)

    def update_enriched_targets(self):  # return queryset of targets between the two primers and updates the m2m targets field
        assert self.left.chromosome == self.right.chromosome
        assert self.chromosome == self.left.chromosome
        self.targets = Target.objects.filter(chromosome=self.chromosome, start_pos__gte=self.left.start_pos)\
            .filter(end_pos__lte=self.right.end_pos)\
            .exclude(pk__in=Primer.objects.all().values('pk'))
        self.save()
        return self.targets.all()

    def amplicon_indices(self):
        return (self.left.start_pos, self.right.end_pos)

    def __unicode__(self):
        return 'TE: left=%s, right=%s' % (self.left.name, self.right.name)
### -------------------------------------------------------------------------------------
class Coordinates(models.Model):#XYZ coordinates of laser capture.
    x = models.DecimalField(max_digits=10, decimal_places=4)
    y = models.DecimalField(max_digits=10, decimal_places=4)
    z = models.DecimalField(max_digits=10, decimal_places=4)
    #TODO: discuss Z vs slides
    def __unicode__(self):
        return '(x={0},y={1},z={2})'.format(self.x,self.y,self.z)
### -------------------------------------------------------------------------------------
class FACSMarker(models.Model):
    name = models.CharField(max_length=50)
### -------------------------------------------------------------------------------------
class Panel(models.Model):#collection of targets
    name = models.CharField(max_length=50)
    targets = models.ManyToManyField(TargetEnrichment, related_name='panels')
    #TODO: add comment field
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Location(models.Model):  # Freetext location
    name = models.CharField(max_length=50)

    def __unicode__(self):#TODO: exapnd
        return self.name
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Algorithms description
### -------------------------------------------------------------------------------------
class AlgorithmType(models.Model):#e.g. raw data to consensus data
    name = models.CharField(max_length=50)
    input = models.CharField(max_length=50)
    output = models.CharField(max_length=50)
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class Algorithm(models.Model):
    name = models.CharField(max_length=50)
    type = models.ForeignKey(AlgorithmType)
    version = models.CharField(max_length=50)
    developers = models.ManyToManyField(User)
    def __unicode__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('algorithm-details', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class AlgorithmParameter(models.Model):
    algorithm = models.ManyToManyField(Algorithm, related_name='parameters')
    name = models.CharField(max_length=50)
    type = models.IntegerField() # todo: choices
    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class AlgorithmRun(models.Model):  # invocation
    algorithm = models.ForeignKey(Algorithm, related_name='runs')
    runname = models.CharField(max_length=50, null=True)
    parameters = models.ManyToManyField(AlgorithmParameter, through='AlgorithmRunParameters')
    user = models.ForeignKey(User)
    timestamp = models.DateTimeField()
    status = models.IntegerField() # todo: choices
    ExtraFiles = models.ManyToManyField(ExperimentFile)
    def __unicode__(self):
        return self.runname
    def experiment(self):
#        if self.algorithm.type == AlgorithmType.objects.get() #TODO: Consider constraining the following 'if's according to AlgorithmType
        if self.crd.all():
            return self.crd.all()[0].experiment()
        if self.targetsdistributions.all():
            return self.targetsdistributions.all()[0].experiment()
        # if self.truesequences.all():
        #     return self.truesequences.all()[0].experiment()
        if self.targetsvariants.all():
            return self.targetsvariants.all()[0].experiment()
        if self.gensigs.all():
            return self.gensigs.all()[0].experiment()
        if self.dms.all():
            return self.dms.all()[0].experiment()
        print 'ERROR: no experiment associated with SequenceDistribution!'
        raise #TODO: implement proper exception here


### -------------------------------------------------------------------------------------
class AlgorithmRunParameters(models.Model):
    run = models.ForeignKey(AlgorithmRun)
    parameter = models.ForeignKey(AlgorithmParameter)
    value = models.CharField(max_length=50)
### -------------------------------------------------------------------------------------
####### TODO: Log #############

### -------------------------------------------------------------------------------------
### Sampling Hierarchy
### -------------------------------------------------------------------------------------
class Individual(models.Model):
    GENDER = (('M', 'Male'), ('F', 'Female'),)
    taxa = models.ForeignKey(Taxa)
    name = models.CharField(max_length=50)
    sex = models.CharField(max_length=1, choices=GENDER)
    born = models.DateTimeField(null=True, blank=True)
    comment = models.TextField(null=True, blank=True)
    background = models.ForeignKey(GeneticBackground, null=True, blank=True)
    location = models.ForeignKey(Location, null=True, blank=True)
    sacrificed = models.DateTimeField(null=True, blank=True)

    def __unicode__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('individual_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class ExtractionEvent(models.Model):
    individual = models.ForeignKey(Individual)
    name = models.CharField(max_length=100)
    comment = models.TextField(null=True, blank=True)
    date = models.DateTimeField()
    location = models.ForeignKey(Location, null=True, blank=True)
    user_performed = models.ForeignKey(User, related_name='+')
    user_documented = models.ForeignKey(User, related_name='+')

    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.individual), self.name)

    def get_absolute_url(self):
        return reverse('extraction_event_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class Extraction(models.Model):
    extraction_event = models.ForeignKey(ExtractionEvent)
    name = models.CharField(max_length=50)
    date = models.DateTimeField(null=True, blank=True)
    organ = models.ForeignKey(Organ)
    tissue = models.ForeignKey(Tissue)
    comment = models.TextField(null=True, blank=True)
    physical_locations = generic.GenericRelation('SampleLocation',
                               content_type_field='content_type',
                               object_id_field='object_id')
    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.extraction_event), self.name)

    def get_absolute_url(self):
        return reverse('extraction_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class SamplingEvent(models.Model):
    def path(self, filename):
        return settings.MEDIA_ROOT+'/CLineageFiles/sampling_events/'+self.name+'/'+filename
    name = models.CharField(max_length=50)
    extraction = models.ForeignKey(Extraction)
    date = models.DateField(null=True, blank=True)
    user = models.ForeignKey(User, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)
    attachment = models.FileField(upload_to=path, null=True, blank=True)
    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.extraction), self.name)

    def get_absolute_url(self):
        return reverse('sampling_event_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class FACS(SamplingEvent):
    marker = models.ForeignKey(FACSMarker)

    def get_absolute_url(self):
        return reverse('facs_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class LaserCapture(SamplingEvent):
    coordinates = models.ForeignKey(Coordinates, null=True, blank=True)  # every cell has distinct coordinates.

    def get_absolute_url(self):
        return reverse('lasercapture_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class CellSelector(SamplingEvent):
    coordinates = models.ForeignKey(Coordinates, null=True, blank=True)  # every cell has distinct coordinates.

    def get_absolute_url(self):
        return reverse('cellselector_detail', kwargs={'pk': self.pk})
### -------------------------------------------------------------------------------------
class Cell(models.Model):
    sampling = models.ForeignKey(SamplingEvent)
    name = models.CharField(max_length=50)
    experiment = models.ManyToManyField(Experiment, related_name='cells', null=True, blank=True)
    composition = models.ForeignKey(SampleComposition)#single cell or bulk
    status = models.ForeignKey(SampleStatus, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)

    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.sampling), self.name)

    def get_absolute_url(self):
        return reverse('cell_detail', kwargs={'pk': self.pk})
## -------------------------------------------------------------------------------------
class CellContent(models.Model):  # aka DNA
    parent = models.ForeignKey('CellContent', null=True, blank=True)
    cell = models.ForeignKey(Cell)
    panel = models.ForeignKey(Panel, null=True, blank=True)
    type = models.ForeignKey(CellContentType)
    name = models.CharField(max_length=50, null=True, blank=True)
    protocol = models.ForeignKey(Protocol, null=True, blank=True)
    seq_ready = models.BooleanField(default=False)
    user = models.ForeignKey(User, null=True, blank=True)
    comment = models.TextField()
    physical_locations = generic.GenericRelation('SampleLocation',
                               content_type_field='content_type',
                               object_id_field='object_id')
    def __unicode__(self):
        return u'{}>{}'.format(unicode(self.cell), self.name)

    def get_absolute_url(self):
        return reverse('cell_content_detail', kwargs={'pk': self.pk})

    def autoname(self):
        if self.seq_ready:
            return "%s_%s_%s" % self.cell.name, self.protocol.initials, 'seqready'
        return "%s_%s" % self.cell.name, self.protocol.initials

    def experiment(self):
        print self.cell.experiment.values('id').annotate(experiment_count=Count('id')).order_by('-experiment_count')
        return Experiment.objects.get(id = self.cell.experiment.values('id').annotate(experiment_count=Count('id')).order_by('-experiment_count')[0]['id'])
### -------------------------------------------------------------------------------------


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
class Sequencing(models.Model):
    samples = models.ManyToManyField(CellContent)
    data = models.ForeignKey('RawData', related_name='sequencing_event', null=True, blank=True)
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


### -------------------------------------------------------------------------------------
### Sequencing Data classes
### -------------------------------------------------------------------------------------
class RawData(models.Model):
    sequencing = models.ForeignKey(Sequencing)
    file = models.FilePathField(blank=True, null=True)
    user = models.ForeignKey(User)
#    def __init__(self, *args, **kwargs):
#        super(RawData, self).__init__(*args, **kwargs)
#        if self.id: #non-empty instance
#            self.file.path = self.directory()
    def __unicode__(self):
        return '%s_raw'%unicode(self.sequencing)

    def directory(self):
        path = '%s/fastq'%(self.sequencing.directory())
        if not os.path.exists(path):
            os.makedirs(path)
        return path
### -------------------------------------------------------------------------------------
class CorrectedRawData(models.Model):
    rawdata = models.ForeignKey(RawData)
    file = models.FilePathField()
    creation = models.ForeignKey(AlgorithmRun, related_name = 'crd')

    def __unicode__(self):
        return '%s_corrected'%unicode(self.rawdata)

    def experiment(self):
        return self.rawdata.experiment()

    def directory(self):
        path = '%s/corrected'%(self.rawdata.sequencing.directory())
        if not os.path.exists(path):
            os.makedirs(path)
        return path

@receiver(post_save, sender=CorrectedRawData)
def on_save(sender, **kwargs):
    kwargs['instance'].directory()
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Targets Hierarchy
### -------------------------------------------------------------------------------------
class FailedTargetValue(models.Model):
    comment = models.TextField()
### -------------------------------------------------------------------------------------
class SequenceDistribution(models.Model):
    target = models.ManyToManyField(Target, through='TargetAnalysis')
    value = models.ForeignKey(Sequence)
    count = models.IntegerField()
    failed = models.ForeignKey(FailedTargetValue, null=True)

    def __unicode__(self):
        return self.value+'-'+str(self.count)

    def experiment(self):
        experiments = []
        for analysis in self.targetanalysis_set.all():
            experiments.append(analysis.sequencingdata.experiment())
        experiment = list(set(experiments))
        if len(experiment) > 1:
            print 'ERROR: multiple experiments associated with SequenceDistribution!'
            raise #TODO: implement proper exception here
        if experiment:
            return experiment[0]
        print 'ERROR: no experiment associated with SequenceDistribution!'
        raise #TODO: implement proper exception here
### -------------------------------------------------------------------------------------
class TargetAnalysis(models.Model):
    target = models.ForeignKey(Target)
    distribution = models.ForeignKey(SequenceDistribution)
    creation = models.ForeignKey(AlgorithmRun, related_name='targetsdistributions')
    sequencingdata = models.ForeignKey(CorrectedRawData)

    def __unicode__(self):
        return self.target.__unicode__()+'_'+self.distribution.__unicode__()

    def experiment(self):
        return self.sequencingdata.experiment()
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
class TargetVariant(models.Model):
    distribution = models.ManyToManyField(SequenceDistribution)
    creation = models.ForeignKey(AlgorithmRun, related_name='targetsvariants')
    type = models.ForeignKey(TargetVariantType)
    value = models.PositiveIntegerField()
    #TODO:Decide if the same target as parent or different one.

    def __unicode__(self):
        return self.type.__unicode__()+'_'+str(self.value)

    def experiment(self):
        return self.ts.experiment()
### -------------------------------------------------------------------------------------
class GenSig(models.Model):
    variants = models.ManyToManyField(TargetVariant)
    creation = models.ForeignKey(AlgorithmRun, related_name='gensigs')
    value = models.TextField() #vector

    def __unicode__(self):
        return self.value

    def experiment(self):
        experiments = []
        for tvariant in self.variants:
            experiments.append(tvariant.experiment())
        experiment = list(set(experiments))
        if len(experiment) > 1:
            print 'ERROR: multiple experiments associated with GenSig!'
            raise #TODO: implement proper exception here
        if experiment:
            return experiment[0]
        print 'ERROR: no experiment associated with GenSig!'
        raise #TODO: implement proper exception here
### -------------------------------------------------------------------------------------
class DM(models.Model):
    cell1 = models.ForeignKey(GenSig, related_name='+')
    cell2 = models.ForeignKey(GenSig, related_name='+')
    distance = models.IntegerField()
    creation = models.ForeignKey(AlgorithmRun, related_name='dms')

    def __unicode__(self):
        return self.cell1.__unicode__()+'_'+self.cell2.__unicode__()+'='+str(self.distance)

    def experiment(self):
        if self.cell1.experiment() == self.cell2.experiment():
            return self.cell1.experiment()
        else:
            print 'ERROR: multiple experiments associated with DM!'
            raise #TODO: implement proper exception here
### -------------------------------------------------------------------------------------
class CellTreeNode(MPTTModel):
    cell = models.ForeignKey(Cell)
    parent = TreeForeignKey('self', null=True, blank = True, related_name='children')
    distance = models.IntegerField()
    def __unicode__(self):
        return self.cell.name
#    class MPTTMeta:
#        order_insertion_by=['name']
### -------------------------------------------------------------------------------------


### -------------------------------------------------------------------------------------
### Optional storage mapping
### -------------------------------------------------------------------------------------
class StorageType(models.Model):
    name = models.CharField(max_length=100, blank=True)
    temperature = models.DecimalField(null=True, blank=True, max_digits=5, decimal_places=1)
    def __unicode__(self):
        return u"%s (temp:%s\xB0c)" % (self.name, str(self.temperature))
### -------------------------------------------------------------------------------------
class StorageBox(models.Model):
    code = models.AutoField(primary_key=True)
    storage_type = models.ForeignKey(StorageType)
    name = models.CharField(max_length=100, blank=True)
    barcode = models.CharField(max_length=20, blank=True)
    def __unicode__(self):
        return u"%s (type:%s)" % (self.name, self.storage_type.name)
### -------------------------------------------------------------------------------------
class PlateContext(models.Model): #The plate's context in use. e.g. pcr
    description = models.CharField(max_length=30, blank=True)
    def __unicode__(self):
        return u"%s" %self.description
### -------------------------------------------------------------------------------------
class PlatePlastica(models.Model): #The plate's physical form. e.g. deepwell square
    code = models.AutoField(primary_key=True)
    description = models.CharField(max_length=30, blank=True)
    rows = models.IntegerField(default=8)
    columns = models.IntegerField(default=12)

    @property
    def wells(self):
        return self.rows*self.columns

    def __unicode__(self):
        return u"%s" %self.description
### -------------------------------------------------------------------------------------
class PlateType(models.Model):
    code = models.AutoField(primary_key=True)
    friendly = models.CharField(max_length=100)
    context = models.ForeignKey(PlateContext, null=True)
    plastic = models.ForeignKey(PlatePlastica, null=True)
    def __unicode__(self):
        return u"%s" % self.friendly
### -------------------------------------------------------------------------------------
class PlateFullException(Exception):
    """
    Attempt to position reagent when plate is already full
    """
    pass

class Plate(models.Model):
    code = models.AutoField(primary_key=True)
    type = models.ForeignKey(PlateType)
    name = models.CharField(max_length=200, blank=True)
    barcode = models.CharField(max_length=20, blank=True)
    timestamp = models.DateField(null=True, blank=True)
    state = models.CharField(max_length=20, blank=True)
    lastusedwell = models.CharField(max_length=4, default='A1')
    def __unicode__(self):
        return self.name
    def get_absolute_url(self):
        return reverse('plate_detail', kwargs={'pk': self.pk})

    def skip_to_first_free_column(self):  # Can't really remember what this function is for. CADMAD heritage.
        first_free_well_index = str2index(self.lastusedwell)
        if first_free_well_index % 8 != 1:
            self.lastusedwell = index2str(first_free_well_index + 8 - (first_free_well_index-1) % 8)

    def place_reagent(self, reagent, reagent_volume=-1, reagent_concentration=-1):
        print str2index(self.lastusedwell), self.type.plastic.wells
        if str2index(self.lastusedwell) > self.type.plastic.wells:
            raise PlateFullException
        location = SampleLocation.objects.create(plate=self,
                                                 well=self.lastusedwell,
                                                 reagent=reagent,
                                                 volume=reagent_volume,
                                                 concentration=reagent_concentration)
        self.lastusedwell = index2str(str2index(self.lastusedwell)+1)  # increment first free well
        return location
    place_reagent.PlateFullException = PlateFullException
### -------------------------------------------------------------------------------------
class PlateStorage(models.Model):
    storageBox = models.ForeignKey(StorageBox)
    plate = models.ForeignKey(Plate)
    inner_location = models.CharField(max_length=100, blank=True)
    notes = models.CharField(max_length=250, blank=True)
    def __unicode__(self):
        return '%s in %s' % (self.plate.name, self.storageBox.name)
### -------------------------------------------------------------------------------------
class SampleLocation(models.Model):
    plate = models.ForeignKey(Plate)
    well = models.CharField(max_length=3, blank=True)
    content_type = models.ForeignKey(ContentType)
    object_id = models.PositiveIntegerField()
    reagent = generic.GenericForeignKey('content_type', 'object_id')
    volume = models.DecimalField(null=True, max_digits=10, decimal_places=3, blank=True)
    concentration = models.DecimalField(null=True, max_digits=10, decimal_places=5, blank=True)
    timestamp = models.DateTimeField(auto_now=True)
    def __unicode__(self):
        return 'plate %s at %s' % (self.plate.name, self.well)
### -------------------------------------------------------------------------------------
### Primers Multiplex
### -------------------------------------------------------------------------------------
class PrimersMultiplex(models.Model):
    name = models.CharField(max_length=20)
    primers = models.ManyToManyField(TargetEnrichment)
    physical_locations = generic.GenericRelation(SampleLocation,
                               content_type_field='content_type',
                               object_id_field='object_id')
    def __unicode__(self):
        return self.name