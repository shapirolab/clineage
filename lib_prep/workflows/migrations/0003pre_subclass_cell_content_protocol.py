# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

from frogress import bar


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0003_add_fks_after_split'),
        ('workflows', '0002_add_wfc'),
    ]

    operations = [
        migrations.RenameField(
            model_name='cellcontent',
            old_name='protocol',
            new_name='old_protocol',
        ),
        # NOTE: this is here because of bug in django deferred sql.
    ]
