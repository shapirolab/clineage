# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('runs', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DemultiplexedReads',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('fastq1', models.FilePathField(null=True)),
                ('fastq2', models.FilePathField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Demultiplexing',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
        ),
        migrations.CreateModel(
            name='DemultiplexingScheme',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='MergedReads',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('assembled_fastq', models.FilePathField(null=True)),
                ('discarded_fastq', models.FilePathField(null=True)),
                ('unassembled_forward_fastq', models.FilePathField(null=True)),
                ('unassembled_reverse_fastq', models.FilePathField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='MergingScheme',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='NGSKit',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
        ),
        migrations.AddField(
            model_name='machinetype',
            name='read_length',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='ngsrun',
            name='bcl_directory',
            field=models.FilePathField(null=True, allow_files=False, allow_folders=True),
        ),
    ]
