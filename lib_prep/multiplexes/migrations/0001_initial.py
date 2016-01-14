# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('planning', '0001_initial'),
        ('reagents', '0001_initial'),
        ('linapp', '0003_add_fks_after_split'),
    ]

    operations = [
        migrations.CreateModel(
            name='Panel',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='PCR1Multiplex',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=20)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='pcr1multiplex',
            name='ters',
            field=models.ManyToManyField(to='reagents.PCR1PrimerPairTERBase'),
        ),
        migrations.AddField(
            model_name='panel',
            name='tes',
            field=models.ManyToManyField(related_name='panels', to='planning.TargetEnrichment'),
        ),
    ]
