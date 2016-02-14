# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0004_split_targets_etc'),
        ('planning', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='ugsplus',
            name='old_primer',
        ),
        migrations.RemoveField(
            model_name='ugsminus',
            name='old_primer',
        ),
        migrations.RemoveField(
            model_name='target',
            name='old_target',
        ),
        migrations.RemoveField(
            model_name='targetenrichment',
            name='old_te',
        ),
    ]
