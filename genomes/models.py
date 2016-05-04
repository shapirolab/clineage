import os
import mmap

from django.db import models
from django.contrib.auth.models import User
from django.conf import settings

#from targeted_enrichment.planning.models import RestrictionSite

from utils.SequenceManipulations import *
from misc.models import Taxa
from misc.dna import DNA

class SearchMarginesDoesNotExist(Exception):
    """
    Attempt to query a position outside chromosome sequence
    """
    pass


class Assembly(models.Model):
    taxa = models.ForeignKey(Taxa)
    name = models.CharField(max_length=50)
    friendly_name = models.CharField(max_length=50)

    def __str__(self):
        return self.friendly_name

    class Meta:
        verbose_name_plural = 'Assemblies'

    def get_path(self):
        return os.path.join(self.taxa.friendly_name, self.friendly_name)


class Chromosome(models.Model):
    name = models.CharField(max_length=50)
    assembly = models.ForeignKey(Assembly)
    sequence_length = models.IntegerField(null=True)
    cyclic = models.BooleanField()

    def __str__(self):
        return "{}:{}".format(self.assembly, self.name)

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
                return mm[start-1:stop].decode("ASCII").upper()
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


class DNASlice(models.Model):
    chromosome = models.ForeignKey(Chromosome)
    # TODO: decide indexing method
    start_pos = models.IntegerField(db_index=True)
    end_pos = models.IntegerField(db_index=True)
    _sequence = models.CharField(max_length=300,null=True,default=None)

    @property
    def sequence(self):
        if self._sequence is not None:
            return DNA(self._sequence)
        else:
            return DNA(self._get_seq())

    def __len__(self):
        return self.end_pos-self.start_pos+1

    def _get_seq(self):
        return self.chromosome.getdna(self.start_pos,self.end_pos)

    def cache(self):
        self._sequence = self._get_seq()
        self.save()

    def __str__(self):
        return "{}:{}-{}".format(self.chromosome.name, self.start_pos,
            self.end_pos)

    def pretty(self, width=2):
        """
        Return a pretty representaion of this slice, with indices and aligned
        prints.
        Width is the length of each line, in blocks of 10.
        """
        ret = ""
        full_width = width*10
        seq = ((" "*((self.start_pos-1)%full_width)) + "{}" \
            + (" "*((width-1)-((self.end_pos-1)%full_width)))).format(
                self.sequence)
        sa = (self.start_pos-1)//full_width
        ea = (self.end_pos-1)//full_width
        ret += ("{pos!s:<12}: "+" ".join(["{ind}"]*width)+"\n").format(
            pos=self.chromosome, ind="1234567890")
        for i in range(ea-sa+1):
            ret += ("{pos:<12}: "+" ".join(["{}"]*width)+"\n").format(
                *[seq[i*full_width+j*10:i*full_width+(j+1)*10] for \
                    j in range(width)],
                pos=(sa+i)*full_width)
        return ret
