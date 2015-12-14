# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('lib_prep', '0002_auto_20151126_1217'),
    ]

    state_ops = [
        migrations.RemoveField(
            model_name='machine',
            name='type',
        ),
        migrations.RemoveField(
            model_name='panel',
            name='targets',
        ),
        migrations.RemoveField(
            model_name='primersmultiplex',
            name='primers',
        ),
        migrations.RemoveField(
            model_name='sequencing',
            name='machine',
        ),
        migrations.RemoveField(
            model_name='sequencing',
            name='user',
        ),
        migrations.DeleteModel(
            name='Machine',
        ),
        migrations.DeleteModel(
            name='MachineType',
        ),
        migrations.DeleteModel(
            name='Panel',
        ),
        migrations.DeleteModel(
            name='PrimersMultiplex',
        ),
        migrations.DeleteModel(
            name='Sequencing',
        ),
    ]

    db_ops = [
        migrations.AlterModelTable(
            name='Machine',
            table='runs_machine',
        ),
        migrations.AlterModelTable(
            name='MachineType',
            table='runs_machinetype',
        ),
        migrations.AlterModelTable(
            name='Sequencing',
            table='runs_ngsrun',
        ),
    ]

    operations = [
        migrations.RemoveField(
            model_name='sequencing',
            name='protocol',
        ),
        migrations.RemoveField(
            model_name='sequencing',
            name='samples',
        ),
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
            database_operations=db_ops,
        ),
    ]
