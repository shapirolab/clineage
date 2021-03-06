# -*- coding: utf-8 -*-
# Generated by Django 1.9.4 on 2016-08-02 18:05
from __future__ import unicode_literals

from Bio import SeqIO
from frogress import bar

from django.db import migrations, models


def populate_num_reads(apps, schema_editor):
    print()
    db_alias = schema_editor.connection.alias
    SampleReads = apps.get_model("analysis", "SampleReads")
    for sr in bar(SampleReads.objects.using(db_alias).all()):
        with open(sr.fastq1) as f:
            seq = SeqIO.parse(f, format="fastq")
            sr.num_reads = sum(1 for x in seq)
        sr.save()


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0007_add_separation_finished_fields'),
    ]

    operations = [
        migrations.AddField(
            model_name='samplereads',
            name='num_reads',
            field=models.PositiveIntegerField(default=0),
            preserve_default=False,
        ),
        migrations.RunPython(
            code=populate_num_reads,
            reverse_code=migrations.RunPython.noop,
        ),
    ]
