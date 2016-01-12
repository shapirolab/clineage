# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

def split_targets_etc(apps, schema_editor):
    split_targets(apps, schema_editor)
    split_target_enrichments(apps, schema_editor)
    split_multiplexes(apps, schema_editor)
    split_panels(apps, schema_editor)

def rejoin_targets_etc(apps, schema_editor):
    rejoin_targets(apps, schema_editor)
    rejoin_restriction_sites(apps, schema_editor)
    rejoin_primers(apps, schema_editor)
    rejoin_target_enrichments(apps, schema_editor)
    rejoin_multiplexes(apps, schema_editor)
    rejoin_panels(apps, schema_editor)

class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0003_add_fks_after_split'),
        ('genomes', '0002_create_dnaslice'),
        ('reagents', '0001_initial'),
        ('synthesis', '0001_initial'),
        ('parts', '0001_initial'),
        ('planning', '0001_initial'),
        ('multiplexes', '0001_initial'),
    ]

    operations = [
        # All of these fields are temporary, to assist with the data migration.
        migrations.AddField(
            model_name='target',
            name='new_target',
            field=models.ForeignKey(to='planning.Target', null=True),
        ),
        # GFK style.
        migrations.AddField(
            model_name='primer',
            name='new_primer_id',
            field=models.PositiveIntegerField(db_index=True,null=True),
        ),
        migrations.AddField(
            model_name='primer',
            name='new_primer_model',
            #name='content_type',
            field=models.CharField(max_length=30,null=True),
            #field=models.ForeignKey(to='contenttypes.ContentType',null=True),
        ),
        migrations.AddField(
            model_name='primer',
            name='new_ugs_id',
            field=models.PositiveIntegerField(db_index=True,null=True),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='new_te',
            field=models.ForeignKey(to='planning.TargetEnrichment', null=True),
        ),
        # GFK style.
        migrations.AddField(
            model_name='targetenrichment',
            name='new_ter_id',
            field=models.PositiveIntegerField(db_index=True,null=True),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='new_ter_model',
            #name='content_type',
            field=models.CharField(max_length=50,null=True),
            #field=models.ForeignKey(to='contenttypes.ContentType',null=True),
        ),
        migrations.RunPython(
            code=split_targets_etc,
            reverse_code=rejoin_targets_etc,
        ),
    ]
