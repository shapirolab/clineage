import os
import mmap

from django.db import models
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


class DNASliceBase(models.Model):
    chromosome = models.ForeignKey(Chromosome)
    # TODO: decide indexing method
    start_pos = models.IntegerField(db_index=True)
    end_pos = models.IntegerField(db_index=True)
    _sequence = models.CharField(max_length=300, null=True, default=None)

    class Meta:
        abstract = True

    @property
    def sequence(self):
        if self._sequence is not None:
            return DNA(self._sequence)
        else:
            return DNA(self._get_seq())

    def __len__(self):
        return self.end_pos - self.start_pos + 1

    def __lt__(self, other):
        if self.chromosome_id != other.chromosome_id:
            return NotImplemented
        if self.start_pos != other.start_pos:
            return self.start_pos < other.start_pos
        return self.end_pos < other.end_pos

    def __gt__(self, other):
        if self.chromosome_id != other.chromosome_id:
            return NotImplemented
        if self.start_pos != other.start_pos:
            return self.start_pos > other.start_pos
        return self.end_pos > other.end_pos

    def _get_seq(self):
        return self.chromosome.getdna(self.start_pos, self.end_pos)

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
        full_width = width * 10
        seq = ((" " * ((self.start_pos - 1) % full_width)) + "{}" \
               + (" " * ((width - 1) - ((self.end_pos - 1) % full_width)))).format(
            self.sequence)
        sa = (self.start_pos - 1) // full_width
        ea = (self.end_pos - 1) // full_width
        ret += ("{pos!s:<12}: " + " ".join(["{ind}"] * width) + "\n").format(
            pos=self.chromosome, ind="1234567890")
        for i in range(ea - sa + 1):
            ret += ("{pos:<12}: " + " ".join(["{}"] * width) + "\n").format(
                *[seq[i * full_width + j * 10:i * full_width + (j + 1) * 10] for \
                  j in range(width)],
                pos=(sa + i) * full_width)
        return ret


class DNASlice(DNASliceBase):
    contains = models.ManyToManyField('genomes.DNASlice',
                                      related_name='contained',
                                      symmetrical=False,
                                      through='genomes.DNASlice_Contains',
                                      through_fields=('outer', 'inner')
                                      )

    overlaps = models.ManyToManyField('genomes.DNASlice',
                                      related_name='+',
                                      symmetrical=False,
                                      through='genomes.DNASlice_Overlaps',
                                      through_fields=('slice1', 'slice2')
                                      )

    def relative_pos(self, chromosome, position):
        """
        The function checks if nucleotide is in the slice and gives it's relative position in the slice
        :param absolute_pos: tuple of chromosome object and start position of the slice
        :return: relative position
        """
        if self.chromosome != chromosome:
            raise ValueError("The nucleotide is not in same chromosome")
        if position < self.start_pos or position > self.end_pos:
            raise ValueError("The nucleotide is not in DNA slice")

        relative_pos = position - self.start_pos
        return relative_pos

    class Meta:
        unique_together = [
            ("chromosome", "start_pos", "end_pos"),
        ]
        index_together = [
            ("chromosome", "start_pos", "end_pos"),
        ]
        ordering = ["chromosome", "start_pos", "end_pos"]


class DNASlice_Contains(models.Model):
    inner = models.ForeignKey(DNASlice, on_delete=models.DO_NOTHING, related_name='+')
    outer = models.ForeignKey(DNASlice, on_delete=models.DO_NOTHING, related_name='+')

    class Meta:
        managed = False


class DNASlice_Overlaps(models.Model):
    slice1 = models.ForeignKey(DNASlice, on_delete=models.DO_NOTHING, related_name='+')
    slice2 = models.ForeignKey(DNASlice, on_delete=models.DO_NOTHING, related_name='+')

    class Meta:
        managed = False


class RestrictionSiteDNASlice(DNASliceBase):
    contains = models.ManyToManyField('genomes.RestrictionSiteDNASlice',
                                      related_name='contained',
                                      symmetrical=False,
                                      through='genomes.RestrictionSiteDNASlice_Contains',
                                      through_fields=('outer', 'inner')
                                      )


    overlaps = models.ManyToManyField('genomes.RestrictionSiteDNASlice',
                                  related_name='+',
                                  symmetrical=False,
                                  through='genomes.RestrictionSiteDNASlice_Overlaps',
                                  through_fields=('slice1', 'slice2')
                                  )

    class Meta:
        unique_together = [
            ("chromosome", "start_pos", "end_pos"),
        ]
        index_together = [
            ("chromosome", "start_pos", "end_pos"),
        ]
        ordering = ["chromosome", "start_pos", "end_pos"]


class RestrictionSiteDNASlice_Contains(models.Model):
    inner = models.ForeignKey(RestrictionSiteDNASlice, on_delete=models.DO_NOTHING, related_name='+')
    outer = models.ForeignKey(RestrictionSiteDNASlice, on_delete=models.DO_NOTHING, related_name='+')

    class Meta:
        managed = False


class RestrictionSiteDNASlice_Overlaps(models.Model):
    slice1 = models.ForeignKey(RestrictionSiteDNASlice, on_delete=models.DO_NOTHING, related_name='+')
    slice2 = models.ForeignKey(RestrictionSiteDNASlice, on_delete=models.DO_NOTHING, related_name='+')

    class Meta:
        managed = False


class DupSlice(models.Model):
    slice = models.ForeignKey(DNASlice)
    strand = models.CharField(max_length=1)
    uid = models.IntegerField()

    def __str__(self):
        return "{} <==> ???".format(self.slice)

    class Meta:
        unique_together=[
            ("uid", "slice")
        ]
        ordering=["uid", "slice"]