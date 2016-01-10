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
        # These are temp fields.
        migrations.RemoveField(
            model_name='targetenrichment',
            name='new_te',
        ),
        migrations.RemoveField(
            model_name='target',
            name='new_target',
        ),
        # Up to here.
        #migrations.RemoveField(
            #model_name='microsatellite',
            #name='target_ptr',
        #),
        migrations.RemoveField(
            model_name='panel',
            name='targets',
        ),
        migrations.RemoveField(
            model_name='primer',
            name='sequence',
        ),
        migrations.RemoveField(
            model_name='primer',
            name='tail',
        ),
        #migrations.RemoveField(
            #model_name='primer',
            #name='target_ptr',
        #),
        migrations.RemoveField(
            model_name='primersmultiplex',
            name='primers',
        ),
        migrations.RemoveField(
            model_name='restrictionsite',
            name='restriction_type',
        ),
        #migrations.RemoveField(
            #model_name='restrictionsite',
            #name='target_ptr',
        #),
        #migrations.RemoveField(
            #model_name='snp',
            #name='target_ptr',
        #),
        migrations.RemoveField(
            model_name='target',
            name='chromosome',
        ),
        migrations.RemoveField(
            model_name='target',
            name='partner',
        ),
        migrations.RemoveField(
            model_name='target',
            name='referencevalue',
        ),
        migrations.RemoveField(
            model_name='target',
            name='type',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='chromosome',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='left',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='partner',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='right',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='targets',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='type',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='validation_failure',
        ),
        migrations.RemoveField(
            model_name='targetenrichmenttype',
            name='protocol',
        ),
        migrations.DeleteModel(
            name='Microsatellite',
        ),
        migrations.DeleteModel(
            name='Panel',
        ),
        migrations.DeleteModel(
            name='Primer',
        ),
        migrations.DeleteModel(
            name='PrimersMultiplex',
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
            name='Sequence',
        ),
        migrations.DeleteModel(
            name='SNP',
        ),
        migrations.DeleteModel(
            name='Target',
        ),
        migrations.DeleteModel(
            name='TargetEnrichment',
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
    ]
