# -*- coding: utf-8 -*-
# Generated manually, in form of Django 1.9.4 on 2016-05-02 19:56
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


def populate_ter_amplicons(apps, schema_editor):
    pass


class Migration(migrations.Migration):

    dependencies = [
        ('amplicons', '0001_initial'),
        ('reagents', '0004_auto_20160420_1723'),
    ]

    operations = [
        migrations.RunPython(
            code=populate_ter_amplicons,
            reverse_code=migrations.RunPython.noop,
        ),
    ]
