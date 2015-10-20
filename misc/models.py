import hashlib
import mmap
import os

from django.db import models


### -------------------------------------------------------------------------------------
### Generic Biological Objects
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

from django.db import models

import utils.SequenceManipulations as sm

class DNA(object):

    def __init__(self,seq):
        self._seq = seq.upper()

    @property
    def seq(self):
        return self._seq

    def rev_comp(self):
        return DNA(sm.complement(self.seq)[::-1])

    def rev(self):
        return DNA(self.seq[::-1])

    def comp(self):
        return DNA(sm.complement(self.seq))

    def __iadd__(self, other):
        if not isinstance(other,DNA):
            raise TypeError("DNA")
        self._seq += other.seq
        return self

    def __getitem__(self,num_or_slice):
        if not isinstance(num_or_slice,slice):
            raise ValueError("Only returns slices")
        if num_or_slice.step is not None:
            raise ValueError("No step")
        return DNA(self.seq[num_or_slice])

    def __str__(self):
        return self.seq
