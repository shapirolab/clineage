# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('multiplexes', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='PCR1MultiplexCollection',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('mpxs', models.ManyToManyField(to='multiplexes.PCR1Multiplex')),
                ('panel', models.ForeignKey(to='multiplexes.Panel')),
            ],
        ),
    ]
