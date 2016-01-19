# -*- coding: utf-8 -*-
import re
import hashlib

from linapp.migrations.mig_0004.utils import getdna

def get_or_create_sequence(seq, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    Sequence = apps.get_model("linapp", "Sequence")
    if not re.match('^[ACTGactg]+$', seq.strip()):
        raise ValueError('unsupported characters in input sequence {}'.format(seq))
    hash = hashlib.md5(seq).hexdigest()
    try:
        sequence = Sequence.objects.using(db_alias) \
            .get(hash=hash)
    except Sequence.DoesNotExist:
        sequence, created = Sequence.objects.using(db_alias) \
            .get_or_create(
                length=len(seq),
                sequence=seq,
                hash=hash,
        )
    return sequence

class FuncCacheDict(dict):
    def __init__(self,f):
        super(dict,self).__init__()
        self._f = f

    def __missing__(self,x):
        ret = self[x] = self._f(*x)
        return ret

def _get_or_create_target_type(name, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    OldTargetType = apps.get_model("linapp", "TargetType")
    type, c = OldTargetType.objects.using(db_alias).get_or_create(name=name)
    return type

target_types = FuncCacheDict(_get_or_create_target_type)

def unpack_slice(slice, apps, schema_editor):
    """
    Create a dict from a slice. Note: requires chromosome, select_related
    when possible.
    """
    # NOTE: we use id because chromosome isn't populated yet.
    d = {k: slice.__dict__[k] for k in [
            "chromosome_id",
            "start_pos",
            "end_pos",
    ]}
    d["referencevalue"] = get_or_create_sequence(
        getdna(slice.chromosome, slice.start_pos, slice.end_pos),
        apps,
        schema_editor,
    )
    return d
