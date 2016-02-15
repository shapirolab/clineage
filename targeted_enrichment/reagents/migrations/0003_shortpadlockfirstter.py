# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('planning', '0003_auto_20160215_1652'),
        ('synthesis', '0003_shortpadlockfirst'),
        ('reagents', '0002_remove_temp_fks'),
    ]

    operations = [
        migrations.CreateModel(
            name='ShortPadlockFirstTER',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(null=True, blank=True)),
                ('comment', models.CharField(max_length=50, null=True, blank=True)),
                ('padlock', models.ForeignKey(to='synthesis.ShortPadlockFirst')),
                ('te', models.ForeignKey(to='planning.TargetEnrichment')),
                ('validation_failure', models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True)),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
