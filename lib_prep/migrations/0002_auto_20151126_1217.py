# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('genomes', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('sampling', '0001_initial'),
        ('lib_prep', '0001_initial'),
    ]

    state_ops = [
        migrations.AddField(
            model_name='sequencing',
            name='samples',
            field=models.ManyToManyField(to='sampling.CellContent'),
        ),
        migrations.AddField(
            model_name='sequencing',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='primersmultiplex',
            name='primers',
            field=models.ManyToManyField(to='genomes.TargetEnrichment'),
        ),
        migrations.AddField(
            model_name='panel',
            name='targets',
            field=models.ManyToManyField(related_name='panels', to='genomes.TargetEnrichment'),
        ),
        migrations.AddField(
            model_name='machine',
            name='type',
            field=models.ForeignKey(to='lib_prep.MachineType'),
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops
        ),
    ]
