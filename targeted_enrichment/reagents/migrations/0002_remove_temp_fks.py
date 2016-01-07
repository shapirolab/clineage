# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0004_split_targets_etc'),
        ('reagents', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='pcr1primerpairterbase',
            name='old_te',
        ),
        migrations.RemoveField(
            model_name='targetednotailprimerpairter',
            name='old_te',
        ),
    ]
