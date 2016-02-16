# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reagents', '0002_remove_temp_fks'),
        ('planning', '0002_remove_temp_fks'),
        ('synthesis', '0002_remove_temp_fks'),
        ('parts', '0002_remove_temp_fks'),
        ('linapp', '0004_split_targets_etc'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='panel',
            name='targets',
        ),
        migrations.RemoveField(
            model_name='primersmultiplex',
            name='primers',
        ),
        migrations.RemoveField(
            model_name='target',
            name='partner',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='targets',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='partner',
        ),
        migrations.DeleteModel(
            name='Microsatellite',
        ),
        migrations.DeleteModel(
            name='Panel',
        ),
        migrations.DeleteModel(
            name='PrimersMultiplex',
        ),
        migrations.DeleteModel(
            name='TargetEnrichment',
        ),
        migrations.DeleteModel(
            name='Primer',
        ),
        migrations.DeleteModel(
            name='PrimerTail',
        ),
        migrations.DeleteModel(
            name='RestrictionSite',
        ),
        migrations.DeleteModel(
            name='RestrictionSiteType',
        ),
        migrations.DeleteModel(
            name='SNP',
        ),
        migrations.DeleteModel(
            name='Target',
        ),
        migrations.DeleteModel(
            name='TargetEnrichmentFailureType',
        ),
        migrations.DeleteModel(
            name='TargetEnrichmentType',
        ),
        migrations.DeleteModel(
            name='TargetType',
        ),
        migrations.DeleteModel(
            name='Sequence',
        ),
    ]
