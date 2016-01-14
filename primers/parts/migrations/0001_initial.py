# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import primers.strand


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0003_add_fks_after_split'),
    ]

    operations = [
        migrations.CreateModel(
            name='DNABarcode1',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='DNABarcode2',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='IlluminaFlowCellAdaptor1',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='IlluminaFlowCellAdaptor2',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='IlluminaReadingAdaptor1',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='IlluminaReadingAdaptor1Cuts',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('overlap_start', models.IntegerField()),
                ('overlap_end', models.IntegerField()),
                ('ira', models.ForeignKey(to='parts.IlluminaReadingAdaptor1')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='IlluminaReadingAdaptor2',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='IlluminaReadingAdaptor2Cuts',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('overlap_start', models.IntegerField()),
                ('overlap_end', models.IntegerField()),
                ('ira', models.ForeignKey(to='parts.IlluminaReadingAdaptor2')),
            ],
            options={
                'abstract': False,
            },
        ),
        # Temporary fields to assist with the data migration.
        migrations.AddField(
            model_name='illuminareadingadaptor1cuts',
            name='old_tail',
            field=models.ForeignKey(to='linapp.PrimerTail', null=True),
        ),
        migrations.AddField(
            model_name='illuminareadingadaptor2cuts',
            name='old_tail',
            field=models.ForeignKey(to='linapp.PrimerTail', null=True),
        ),
    ]
