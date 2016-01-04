# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('synthesis', '0001_initial'),
        ('planning', '0001_initial'),
        ('linapp', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='PCR1PrimerPairTERBase',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(null=True, blank=True)),
                ('comment', models.CharField(max_length=50, null=True, blank=True)),
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
        migrations.CreateModel(
            name='PCR1PrimerPairTER',
            fields=[
                ('pcr1primerpairterbase_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='reagents.PCR1PrimerPairTERBase')),
                ('left_primer', models.ForeignKey(to='synthesis.PCR1PlusPrimer')),
                ('right_primer', models.ForeignKey(to='synthesis.PCR1MinusPrimer')),
            ],
            options={
                'abstract': False,
            },
            bases=('reagents.pcr1primerpairterbase',),
        ),
        migrations.CreateModel(
            name='PCR1PrimerPairTERDeprecated',
            fields=[
                ('pcr1primerpairterbase_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='reagents.PCR1PrimerPairTERBase')),
                ('left_primer', models.ForeignKey(to='synthesis.PCR1PlusPrimer')),
                ('right_primer', models.ForeignKey(to='synthesis.PCR1MinusPrimer')),
            ],
            options={
                'abstract': False,
            },
            bases=('reagents.pcr1primerpairterbase',),
        ),
        migrations.CreateModel(
            name='PCR1WithCompanyTagPrimerPairTER',
            fields=[
                ('pcr1primerpairterbase_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='reagents.PCR1PrimerPairTERBase')),
                ('left_primer', models.ForeignKey(to='synthesis.PCR1WithCompanyTagPlusPrimer')),
                ('right_primer', models.ForeignKey(to='synthesis.PCR1WithCompanyTagMinusPrimer')),
            ],
            options={
                'abstract': False,
            },
            bases=('reagents.pcr1primerpairterbase',),
        ),
        migrations.AddField(
            model_name='targetednotailprimerpairter',
            name='validation_failure',
            field=models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True),
        ),
        migrations.AddField(
            model_name='pcr1primerpairterbase',
            name='te',
            field=models.ForeignKey(to='planning.TargetEnrichment'),
        ),
        migrations.AddField(
            model_name='pcr1primerpairterbase',
            name='validation_failure',
            field=models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True),
        ),
        # Temporary fields to assist with the data migration.
        migrations.AddField(
            model_name='targetednotailprimerpairter',
            name='old_te',
            field=models.ForeignKey(to='linapp.TargetEnrichment'),
        ),
        migrations.AddField(
            model_name='pcr1primerpairterbase',
            name='old_te',
            field=models.ForeignKey(to='linapp.TargetEnrichment'),
        ),
    ]
