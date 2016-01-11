# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

def split_targets_etc(apps, schema_editor):
    print
    split_targets(apps, schema_editor)
    split_target_enrichments(apps, schema_editor)
    split_multiplexes(apps, schema_editor)
    split_panels(apps, schema_editor)

def rejoin_targets_etc(apps, schema_editor):
    print
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
        # All of these fields, up to the RunPython, are temporary, to assist
        # with the data migration.
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
        migrations.RemoveField(
            model_name='microsatellite',
            name='target_ptr',
        ),
        migrations.RemoveField(
            model_name='primer',
            name='sequence',
        ),
        migrations.RemoveField(
            model_name='primer',
            name='tail',
        ),
        migrations.RemoveField(
            model_name='primer',
            name='target_ptr',
        ),
        migrations.RemoveField(
            model_name='restrictionsite',
            name='restriction_type',
        ),
        migrations.RemoveField(
            model_name='restrictionsite',
            name='target_ptr',
        ),
        migrations.RemoveField(
            model_name='snp',
            name='target_ptr',
        ),
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
