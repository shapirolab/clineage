import os
import mmap

from django.db import models
from django.contrib.contenttypes import fields
from django.contrib.auth.models import User
from django.conf import settings

from utils.SequenceManipulations import *
from misc.models import Taxa
from linapp.models import Protocol
from wet_storage.models import SampleLocation

class SearchMarginesDoesNotExist(Exception):
    """
    Attempt to query a position outside chromosome sequence
    """
    pass

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
### Full Genomes
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
        if start > 0 \
           and stop > 0 \
           and start <= stop \
           and stop <= self.sequence_length:
            with open(self.get_abs_path(), 'r+b') as f:
                mm = mmap.mmap(f.fileno(), 0)
                return mm[start-1:stop].upper()
        if self.cyclic:
            if start > self.sequence_length:
                start = start-self.sequence_length
            if stop > self.sequence_length:
                stop = stop-self.sequence_length
            if start <= 0:
                start = self.sequence_length + start
            if stop <= 0:
                stop = self.sequence_length + stop
            if start > stop:
                return self.getdna(start, self.sequence_length) + self.getdna(1, stop)
            return self.getdna(start, stop)

        raise ValueError('indices out of bounds')

    def locate(self, start, stop, sequence, padding=10):
        if self.cyclic:
            l = padding
            r = padding
        else:
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
### Types and descriptors
### -------------------------------------------------------------------------------------
class TargetType(models.Model):
    name = models.CharField(max_length=50)

    def __unicode__(self):
        return self.name
### -------------------------------------------------------------------------------------
class RestrictionSiteType(models.Model):
    name = models.CharField(max_length=50)
    sequence = models.CharField(max_length=50)
    cut_delta = models.IntegerField()  # position of cutting site relative to start_pos
    sticky_bases = models.IntegerField()
    sequence_len = models.PositiveIntegerField()

    def save(self, *args, **kwargs):
        self.sequence_len = len(self.sequence)
        return super(RestrictionSiteType, self).save(*args, **kwargs)

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
### Targets?
### -------------------------------------------------------------------------------------
class Target(models.Model):#Target is a locus on a reference genome.
    name = models.CharField(max_length=50)
    type = models.ForeignKey(TargetType) #Microsatellite / SNP / etc...
    chromosome = models.ForeignKey(Chromosome)
    start_pos = models.IntegerField(db_index=True)
    end_pos = models.IntegerField(db_index=True)
    referencevalue = models.ForeignKey(Sequence)
    partner = models.ManyToManyField(User, null=True)

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

    def get_margine(self, pos):
        if self.chromosome.cyclic:
            return pos
        if 0 <= pos <= self.chromosome.sequence_length:
            return pos
        raise SearchMarginesDoesNotExist


    def get_left_surrounding_restriction(self, restriction_type, max_seek=100):
        left_restriction_site = RestrictionSite.objects.filter(restriction_type=restriction_type)\
            .filter(chromosome=self.chromosome).\
            filter(start_pos__lte=self.start_pos-restriction_type.cut_delta).\
            filter(start_pos__gte=self.get_margine(self.start_pos-max_seek)).order_by('-start_pos')
        if left_restriction_site:
            return left_restriction_site[0]
        return self.get_left_surrounding_restriction(restriction_type, max_seek=max_seek*2)

    def get_right_surrounding_restriction(self, restriction_type, max_seek=100):
        right_restriction_site = RestrictionSite.objects.filter(restriction_type=restriction_type)\
            .filter(chromosome=self.chromosome).\
            filter(end_pos__gte=self.end_pos+(restriction_type.sequence_len-restriction_type.cut_delta)).\
            filter(end_pos__lte=self.get_margine(self.end_pos+max_seek)).order_by('start_pos')
        if right_restriction_site:
            return right_restriction_site[0]
        return self.get_right_surrounding_restriction(restriction_type, max_seek=max_seek*2)

    def get_surrounding_restriction(self, restriction_type):
        left = self.get_left_surrounding_restriction(restriction_type)
        right = self.get_right_surrounding_restriction(restriction_type)
        return left, right
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
    physical_locations = fields.GenericRelation('SampleLocation',
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
### -------------------------------------------------------------------------------------
class SNP(Target):
    mutation = models.CharField(max_length=10) #X>Y
    modified = models.CharField(max_length=10) #Y
### -------------------------------------------------------------------------------------
class RestrictionSite(Target):
    restriction_type = models.ForeignKey(RestrictionSiteType)
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
    physical_locations = fields.GenericRelation('SampleLocation',
                                                 content_type_field='content_type',
                                                 object_id_field='object_id')
    targets = models.ManyToManyField(Target, related_name='primer_pair', null=True, blank=True)
    partner = models.ManyToManyField(User, null=True)

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

    def get_internal_restriction(self, restriction):
        return [self.amplicon_indices()[0] + m.start() for m in re.finditer(restriction, self.chromosome.getdna(*self.amplicon_indices()))]

    def get_surrounding_restriction(self, restriction, max_seek=5000):
        for x in range(0, max_seek, 10):
            lamplicon = self.chromosome.getdna(self.amplicon_indices()[0]-x, self.amplicon_indices()[0])
            lttaas = [self.amplicon_indices()[0] - m.start() for m in re.finditer(restriction, lamplicon)]
            if lttaas:
                break

        for x in range(0, max_seek, 10):
            ramplicon = self.chromosome.getdna(self.amplicon_indices()[1], self.amplicon_indices()[1]+x)
            rttaas = [self.amplicon_indices()[1] + m.start() for m in re.finditer(restriction, ramplicon)]
            if rttaas:
                break

        if lttaas and rttaas:
            return max(lttaas), min(rttaas)
        return None


    def __unicode__(self):
        return 'TE: left=%s, right=%s' % (self.left.name, self.right.name)
### -------------------------------------------------------------------------------------

### -------------------------------------------------------------------------------------
### New
### -------------------------------------------------------------------------------------
class DNABarcode(models.Model):
    #TODO: add boundries for actual DNA content.
    sequencing_primer = models.ForeignKey(SequencingPrimer)
    barcode = models.ForeignKey(Sequence)
    adaptor = models.ForeignKey(Adaptor)
    physical_locations = fields.GenericRelation('SampleLocation',
                                             content_type_field='content_type',
                                             object_id_field='object_id')

### -------------------------------------------------------------------------------------
class SequencingPrimer(models.Model):
    value = models.ForeignKey(Sequence)
### -------------------------------------------------------------------------------------
class Adaptor(models.Model):
    value = models.ForeignKey(Sequence)


