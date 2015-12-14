# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

from frogress import bar

def populate_cell_content_protocol(apps, schema_editor):
    print
    CellContent = apps.get_model("sampling", "CellContent")
    CellContentProtocol = apps.get_model("sampling", "CellContentProtocol")
    db_alias = schema_editor.connection.alias
    for cc in bar(CellContent.objects.using(db_alias).all()):
        protocol = cc.old_protocol
        if protocol is None:
            continue
        try:
            cell_protocol = protocol.cellcontentprotocol
        except CellContentProtocol.DoesNotExist:
            cell_protocol = CellContentProtocol(protocol_ptr_id=protocol.pk)
            cell_protocol.__dict__.update(protocol.__dict__)
            cell_protocol.save(using=db_alias)
        cc.protocol = cell_protocol
        cc.save(using=db_alias)

def populate_protocol(apps, schema_editor):
    print
    CellContent = apps.get_model("sampling", "CellContent")
    db_alias = schema_editor.connection.alias
    for cc in bar(CellContent.objects.using(db_alias).all()):
        protocol = cc.protocol
        if protocol is None:
            continue
        cc.old_protocol = protocol.protocol_ptr
        cc.save(using=db_alias)


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0003_auto_20151127_0133'),
        ('sampling', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='cellcontent',
            name='panel',
        ),
        migrations.RemoveField(
            model_name='cellcontent',
            name='parent',
        ),
        migrations.RemoveField(
            model_name='cellcontent',
            name='seq_ready',
        ),
        migrations.RenameField(
            model_name='cellcontent',
            old_name='protocol',
            new_name='old_protocol',
        ),
    ]
