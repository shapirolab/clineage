# -*- coding: utf-8 -*-
import os

from django.db import transaction

class NotCovering(Exception):
    pass

class NotMutuallyExclusive(Exception):
    pass

def get_physical_locations(old_obj, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    ContentType = apps.get_model("contenttypes", "ContentType")
    SampleLocation = apps.get_model("wet_storage", "SampleLocation")
    # NOTE: this relies on get_for_model creating the new types if required,
    # since we don't have the new CTs available yet.
    old_ct = ContentType.objects.get_for_model(old_obj)
    locations = SampleLocation.objects.using(db_alias).filter(
        content_type=old_ct,
        object_id=old_obj.id,
    )
    return locations

def alter_physical_locations(new_obj, locations, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    ContentType = apps.get_model("contenttypes", "ContentType")
    new_ct = ContentType.objects.get_for_model(new_obj)
    for location in locations:
        timestamp = [f for f in location._meta.fields if f.name=="timestamp"][0]
        timestamp.auto_now = False
        location.content_type = new_ct
        location.object_id = new_obj.id
        location.save()
        timestamp.auto_now = True

def transfer_physical_locations(old_obj, new_obj, apps, schema_editor):
    """
    Transfer all physical locations (GFK in SampleLocation-s) pointing at
    old_obj to point at new_obj.
    """
    locations = get_physical_locations(old_obj, apps, schema_editor)
    alter_physical_locations(new_obj, locations, apps, schema_editor)

def get_partner_ids(obj, apps, schema_editor):
    """
    Assuming obj has an M2M to User called "partner", selects all ids
    related to it, without selecting the users themselves.
    """
    db_alias = schema_editor.connection.alias
    return [oop.user_id for oop in obj.partner.through.objects \
        .using(db_alias).filter(**{obj._meta.model_name: obj})]

def check_partition(related_manager, queries):
    """
    Check that the given list of queries (Q objects or kwargs dicts) partitions
    the given related_manager: they are mutually exclusive, and cover the entire
    manager.
    """
    not_covering = related_manager
    querysets = [None] * len(queries)
    for i,q in enumerate(queries):
        if isinstance(q,dict):
            not_covering = not_covering.exclude(**q)
            querysets[i] = related_manager.filter(**q)
        elif isinstance(q,models.Q):
            not_covering = not_covering.exclude(q)
            querysets[i] = related_manager.filter(q)
        else:
            raise TypeError("Bad query type.")
    if not_covering:
        raise NotCovering(not_covering)
    for i in range(len(querysets)):
        for j in range(i):
            if isinstance(queries[j],dict):
                mutex = querysets[i].filter(**queries[j])
            else: # We already checked types above.
                mutex = querysets[i].filter(queries[j])
            if mutex:
                raise NotMutuallyExclusive(i,j,mutex)

# Copypasta from current Assembly & Chromosome.
from django.conf import settings
import mmap

def _assembly_get_path(self):
    return os.path.join(self.taxa.friendly_name, self.friendly_name)

def _chromosome_get_path(self, ext="txt"):
    return os.path.join(_assembly_get_path(self.assembly), 'chr{}.{}'.format(self.name, ext))

def _chromosome_get_abs_path(self):
    return os.path.join(settings.CHROMOSOMES_PATH, _chromosome_get_path(self))

# 1-based, edge included indeces (fences instead of fenceposts).
def getdna(self, start, stop):
    if start > 0 \
       and stop > 0 \
       and start <= stop \
       and stop <= self.sequence_length:
        with open(_chromosome_get_abs_path(self), 'r+b') as f:
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
            return getdna(self, start, self.sequence_length) + getdna(self, 1, stop)
        return getdna(self, start, stop)

    raise ValueError('indices out of bounds')

# End copypasta

# Copypasta from current utils.SequenceManipulations
import re 

_basecomp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'t', 't':'a', 'c':'g', 'g':'c'}
def _complement(seq):
    assert re.match('^[ACTGactg]+$',seq)
    compseq = ''
    for base in seq:
        compseq += _basecomp[base]
    return compseq

# End copypasta

def rc(seq):
    return _complement(seq)[::-1]
