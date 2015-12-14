# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('planning', '0001_initial'),
        ('synthesis', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='PCR1PrimerPairTER',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(null=True, blank=True)),
                ('comment', models.CharField(max_length=50, null=True, blank=True)),
                ('left_primer', models.ForeignKey(to='synthesis.PCR1PlusPrimer')),
                ('right_primer', models.ForeignKey(to='synthesis.PCR1MinusPrimer')),
                ('te', models.ForeignKey(to='planning.TargetEnrichment')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='PCR1PrimerPairTERDeprecated',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(null=True, blank=True)),
                ('comment', models.CharField(max_length=50, null=True, blank=True)),
                ('left_primer', models.ForeignKey(to='synthesis.PCR1PlusPrimer')),
                ('right_primer', models.ForeignKey(to='synthesis.PCR1MinusPrimer')),
                ('te', models.ForeignKey(to='planning.TargetEnrichment')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='TargetedNoTailPrimerPairTER',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(null=True, blank=True)),
                ('comment', models.CharField(max_length=50, null=True, blank=True)),
                ('left_primer', models.ForeignKey(to='synthesis.TargetedNoTailPlusPrimer')),
                ('right_primer', models.ForeignKey(to='synthesis.TargetedNoTailMinusPrimer')),
                ('te', models.ForeignKey(to='planning.TargetEnrichment')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='TargetEnrichmentFailureType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField(null=True, blank=True)),
            ],
        ),
        migrations.AddField(
            model_name='targetednotailprimerpairter',
            name='validation_failure',
            field=models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True),
        ),
        migrations.AddField(
            model_name='pcr1primerpairterdeprecated',
            name='validation_failure',
            field=models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True),
        ),
        migrations.AddField(
            model_name='pcr1primerpairter',
            name='validation_failure',
            field=models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True),
        ),
    ]
