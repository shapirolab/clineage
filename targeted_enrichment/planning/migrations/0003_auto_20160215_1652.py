# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('planning', '0002_remove_temp_fks'),
    ]

    operations = [
        migrations.RenameField(
            model_name='restrictionenzyme',
            old_name='sequence',
            new_name='_sequence',
        ),
    ]
