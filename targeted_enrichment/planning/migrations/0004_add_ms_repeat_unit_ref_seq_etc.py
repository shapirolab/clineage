# -*- coding: utf-8 -*-
# Generated by Django 1.9.4 on 2016-04-20 15:00
from __future__ import unicode_literals

from misc.dna import DNA
from django.db import migrations, models
import django.db.models.deletion
from django.db import IntegrityError
from utils.ms_utils import ms_type_permutations
from frogress import bar
import os


def get_sequence(s):
    if s._sequence is not None:
        return DNA(s._sequence)
    else:
        return DNA(getdna(s.chromosome, s.start_pos, s.end_pos))


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


def get_repeat_unit_ref_seq_forward(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    Microsatellite = apps.get_model("planning", "Microsatellite")
    rut_canon = {}
    print()
    print("Populating repeat_unit_ref_seq in Microsatellite-s:")
    for ms in bar(Microsatellite.objects.using(db_alias).all()):
        repeat_unit_ref_seq = get_sequence(ms.slice)[:ms.repeat_unit_len]
        if repeat_unit_ref_seq not in rut_canon:
            repeat_type_variants = set(ms_type_permutations(repeat_unit_ref_seq))
            canon = min(repeat_type_variants, key=lambda x: x.seq)
            for in_rut in repeat_type_variants:
                rut_canon[in_rut] = canon
        canon = rut_canon[repeat_unit_ref_seq]
        msrut = DNA(ms.repeat_unit_type)
        if msrut != rut_canon[msrut]:
            raise IntegrityError(
                "MS {} has non canonic repeat_unit_type".format(ms))
        if canon != msrut:
            raise IntegrityError(
                "ref seq {} does not match repeat unit type of {}".format(
                    repeat_unit_ref_seq, ms))
        ms.repeat_unit_ref_seq = repeat_unit_ref_seq
        ms.save()


def populate_te_targets(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    Target = apps.get_model("planning", "Target")
    TargetEnrichment = apps.get_model("planning", "TargetEnrichment")
    print()
    print("Populating target M2M in TEs:")
    for te in bar(TargetEnrichment.objects.using(db_alias).all()):
        if te.left.slice.chromosome != te.right.slice.chromosome:
            raise IntegrityError(
                "TE {} has slices from different chromosomes.".format(te))
        if te.chromosome != te.left.slice.chromosome:
            raise IntegrityError(
                "TE {} has a different chromosome than its slices.".format(te))
        # TODO: change to slice query.
        te.targets = Target.objects.filter(
            slice__chromosome=te.chromosome,
            slice__start_pos__gte=te.left.slice.start_pos,
            slice__end_pos__lte=te.right.slice.end_pos,
        )


class Migration(migrations.Migration):

    dependencies = [
        ('genomes', '0003_dnaslice_index'),
        ('planning', '0003_auto_20160215_1652'),
    ]

    operations = [
        migrations.AddField(
            model_name='microsatellite',
            name='planning_version',
            field=models.IntegerField(default=0),
            preserve_default=False,
        ),
        migrations.RunPython(
            code=migrations.RunPython.noop,
            reverse_code=populate_te_targets,
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='targets',
        ),
        migrations.AddField(
            model_name='microsatellite',
            name='repeat_unit_ref_seq',
            field=models.CharField(default='X', max_length=50),
            preserve_default=False,
        ),
        migrations.RunPython(
            code=get_repeat_unit_ref_seq_forward,
            reverse_code=migrations.RunPython.noop,
        ),
        migrations.AlterField(
            model_name='targetenrichment',
            name='chromosome',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='genomes.Chromosome'),
        ),
    ]
