# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('sampling', '0001_initial'),
        ('genomes', '0001_initial'),
        ('misc', '0001_initial'),
        ('runs', '0001_initial'),
        ('wet_storage', '0001_initial'),
        ('workflows', '0001_initial'),
        ('linapp', '0002_split_to_apps'),
    ]

    state_ops = [
        # Readding fields removed because their targets moved.
        migrations.AddField(
            model_name='target',
            name='chromosome',
            field=models.ForeignKey(to='genomes.Chromosome'),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='chromosome',
            field=models.ForeignKey(to='genomes.Chromosome'),
        ),
        migrations.AddField(
            model_name='userreport',
            name='cells',
            field=models.ManyToManyField(to='sampling.Cell'),
        ),
        migrations.AddField(
            model_name='userreport',
            name='individual',
            field=models.ManyToManyField(to='sampling.Individual'),
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
