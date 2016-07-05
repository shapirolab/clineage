# -*- coding: utf-8 -*-
# Generated by Django 1.9.4 on 2016-07-05 18:51
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('planning', '0004_add_ms_repeat_unit_ref_seq_etc'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ugsminus',
            name='slice',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genomes.DNASlice', unique=True),
        ),
        migrations.AlterField(
            model_name='ugsplus',
            name='slice',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genomes.DNASlice', unique=True),
        ),
    ]
