# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0002_split_to_apps'),
    ]

    state_ops = [
        migrations.CreateModel(
            name='Machine',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('machineid', models.CharField(max_length=50)),
                ('name', models.CharField(max_length=100, null=True, blank=True)),
                ('ip', models.IPAddressField(null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='MachineType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('company', models.CharField(max_length=50)),
                ('model', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='Panel',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='PrimersMultiplex',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=20)),
            ],
        ),
        migrations.CreateModel(
            name='Sequencing',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=100)),
                ('date', models.DateField()),
                ('machine', models.ForeignKey(to='lib_prep.Machine')),
                ('protocol', models.ForeignKey(to='linapp.Protocol')),
            ],
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
