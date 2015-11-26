# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0002_split_to_apps'),
    ]

    state_ops = [
        migrations.CreateModel(
            name='Taxa',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('taxonomy_id', models.IntegerField()),
                ('rank', models.CharField(max_length=50)),
                ('parent', models.IntegerField(null=True, blank=True)),
                ('friendly_name', models.CharField(max_length=50)),
            ],
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
