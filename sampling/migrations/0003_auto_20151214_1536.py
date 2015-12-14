# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genomes', '0003_split_targets_etc'),
        ('lib_prep', '0003_kill'),
        ('sampling', '0002_auto_20151130_1421'),
    ]

    state_ops = [
        migrations.RemoveField(
            model_name='cellcontent',
            name='cell',
        ),
        migrations.RemoveField(
            model_name='cellcontent',
            name='protocol',
        ),
        migrations.RemoveField(
            model_name='cellcontent',
            name='type',
        ),
        migrations.RemoveField(
            model_name='cellcontent',
            name='user',
        ),
        migrations.RemoveField(
            model_name='cellcontentprotocol',
            name='protocol_ptr',
        ),
        migrations.DeleteModel(
            name='CellContent',
        ),
        migrations.DeleteModel(
            name='CellContentProtocol',
        ),
        migrations.DeleteModel(
            name='CellContentType',
        ),
    ]

    db_ops = [
        migrations.AlterModelTable(
            name='CellContent',
            table='workflows_cellcontent',
        ),
        migrations.AlterModelTable(
            name='CellContentProtocol',
            table='workflows_cellcontentprotocol',
        ),
        migrations.AlterModelTable(
            name='CellContentType',
            table='workflows_cellcontenttype',
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
            database_operations=db_ops,
        ),
    ]
