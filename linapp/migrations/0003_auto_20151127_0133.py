# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('sampling', '0001_initial'),
        ('linapp', '0002_split_to_apps'),
    ]

    state_ops = [
        migrations.AddField(
            model_name='userreport',
            name='cells',
            field=models.ManyToManyField(to='sampling.Cell'),
        ),
        migrations.AddField(
            model_name='userreport',
            name='individual',
            field=models.ManyToManyField(to='sampling.Individual', null=True),
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]