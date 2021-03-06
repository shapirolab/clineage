# -*- coding: utf-8 -*-
# Generated by Django 1.9.4 on 2016-05-26 17:35
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='adamampliconreads',
            name='fastq1',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adamampliconreads',
            name='fastq2',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adamampliconreads',
            name='fastqm',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adamhistogram',
            name='assignment_sam',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adammarginassignment',
            name='assignment_sam',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adammergedreads',
            name='assembled_fastq',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adammergedreads',
            name='discarded_fastq',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adammergedreads',
            name='unassembled_forward_fastq',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adammergedreads',
            name='unassembled_reverse_fastq',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='adammsvariations',
            name='index_dump_dir',
            field=models.FilePathField(allow_files=False, allow_folders=True, max_length=200),
        ),
        migrations.AlterField(
            model_name='adamreadsindex',
            name='index_dump_dir',
            field=models.FilePathField(allow_files=False, allow_folders=True, max_length=200),
        ),
        migrations.AlterField(
            model_name='histogramentryreads',
            name='fastq1',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='histogramentryreads',
            name='fastq2',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='histogramentryreads',
            name='fastqm',
            field=models.FilePathField(max_length=200, null=True),
        ),
        migrations.AlterField(
            model_name='samplereads',
            name='fastq1',
            field=models.FilePathField(max_length=200),
        ),
        migrations.AlterField(
            model_name='samplereads',
            name='fastq2',
            field=models.FilePathField(max_length=200),
        ),
    ]
