# -*- coding: utf-8 -*-
# Generated by Django 1.9.8 on 2019-02-12 08:17
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('calling', '0007_highestpeaksproximityratiofilteredbisimcorschememodeldot'),
    ]

    operations = [
        migrations.CreateModel(
            name='HighestPeaksMonoSimCorSchemeModelDot',
            fields=[
                ('highestpeaksmonosimcorschememodel_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='calling.HighestPeaksMonoSimCorSchemeModel')),
            ],
            options={
                'abstract': False,
            },
            bases=('calling.highestpeaksmonosimcorschememodel',),
        ),
    ]
