# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('multiplexes', '0001_initial'),
        ('planning', '0001_initial'),
        ('reagents', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='pcr1multiplex',
            name='primers',
            field=models.ManyToManyField(to='reagents.PCR1PrimerPairTERBase'),
        ),
        migrations.AddField(
            model_name='panel',
            name='targets',
            field=models.ManyToManyField(related_name='panels', to='planning.TargetEnrichment'),
        ),
    ]
