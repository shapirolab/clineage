import os
import mmap

from django.db import models
from django.contrib.auth.models import User
from django.conf import settings

from targeted_enrichment.planning.models import RestrictionSite

from utils.SequenceManipulations import *
from misc.models import Taxa
from linapp.models import Protocol

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

    # TODO: decide indexing method
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

###
# New!
###
class DNASlice(models.Model):
    chromosome = models.ForeignKey(Chromosome)
    # TODO: decide indexing method
    start_pos = models.IntegerField(db_index=True)
    end_pos = models.IntegerField(db_index=True)
    _sequence = models.ForeignKey(Sequence,null=True,default=None)

    @property
    def sequence(self):
        if self._sequence is not None:
            return self._sequence
        else:
            return self._get_seq()

    def _get_seq(self):
        self.chromosome.get_dna(self.start_pos,self.end_pos)

    def cache(self):
        self._sequence = self._get_seq()
        self.save()

    def get_margine(self, pos):
        if self.chromosome.cyclic:
            return pos
        if 0 <= pos <= self.chromosome.sequence_length:
            return pos
        raise SearchMarginesDoesNotExist()

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

