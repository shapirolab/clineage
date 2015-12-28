# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import primers.synthesis.models
import primers.strand


class Migration(migrations.Migration):

    dependencies = [
        ('planning', '0001_initial'),
        ('parts', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='PCR1MinusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('irac', models.ForeignKey(to='parts.IlluminaReadingAdaptor2Cuts')),
                ('ugs', models.ForeignKey(to='planning.UGSMinus')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.TargetedHeadMixin, primers.synthesis.models.PCR1TailMixin, models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PCR1PlusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('irac', models.ForeignKey(to='parts.IlluminaReadingAdaptor1Cuts')),
                ('ugs', models.ForeignKey(to='planning.UGSPlus')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.TargetedHeadMixin, primers.synthesis.models.PCR1TailMixin, models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PCR1WithCompanyTagMinusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('irac', models.ForeignKey(to='parts.IlluminaReadingAdaptor2Cuts')),
                ('tag', models.CharField(max_length=1)),
                ('ugs', models.ForeignKey(to='planning.UGSMinus')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.TargetedHeadMixin, primers.synthesis.models.PCR1WithCompanyTagTailMixin, models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PCR1WithCompanyTagPlusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('irac', models.ForeignKey(to='parts.IlluminaReadingAdaptor1Cuts')),
                ('tag', models.CharField(max_length=1)),
                ('ugs', models.ForeignKey(to='planning.UGSPlus')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.TargetedHeadMixin, primers.synthesis.models.PCR1WithCompanyTagTailMixin, models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PCR2MinusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('barcode', models.ForeignKey(to='parts.DNABarcode2')),
                ('ifca', models.ForeignKey(to='parts.IlluminaFlowCellAdaptor2')),
                ('irac', models.ForeignKey(to='parts.IlluminaReadingAdaptor2Cuts')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.PCR2Mixin, primers.strand.MinusStrandMixin, models.Model, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PCR2PlusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('barcode', models.ForeignKey(to='parts.DNABarcode1')),
                ('ifca', models.ForeignKey(to='parts.IlluminaFlowCellAdaptor1')),
                ('irac', models.ForeignKey(to='parts.IlluminaReadingAdaptor1Cuts')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.PCR2Mixin, primers.strand.PlusStrandMixin, models.Model, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='TargetedNoTailMinusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('ugs', models.ForeignKey(to='planning.UGSMinus')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.TargetedHeadMixin, primers.strand.MinusStrandMixin, models.Model, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='TargetedNoTailPlusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('ugs', models.ForeignKey(to='planning.UGSPlus')),
            ],
            options={
                'abstract': False,
            },
            bases=(primers.synthesis.models.TargetedHeadMixin, primers.strand.PlusStrandMixin, models.Model, primers.strand.BaseStrandMixin),
        ),
        # Temporary fields to assist with the data migration.
        migrations.AddField(
            model_name='pcr1plusprimer',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer'),
        ),
        migrations.AddField(
            model_name='pcr1minusprimer',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer'),
        ),
        migrations.AddField(
            model_name='pcr1withcompanytagplusprimer',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer'),
        ),
        migrations.AddField(
            model_name='pcr1withcompanytagminusprimer',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer'),
        ),
        migrations.AddField(
            model_name='targetednotailplusprimer',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer'),
        ),
        migrations.AddField(
            model_name='targetednotailminusprimer',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer'),
        ),
    ]
