# -*- coding: utf-8 -*-
from frogress import bar

from django.db import transaction

from linapp.migrations.mig_0004.utils import transfer_physical_locations

@transaction.atomic
def rejoin_multiplexes(apps, schema_editor):
    print()
    print("Reverting PCR1 multiplexes:")
    db_alias = schema_editor.connection.alias
    PCR1Multiplex = apps.get_model("multiplexes","PCR1Multiplex")
    OldPrimersMultiplex = apps.get_model("linapp","PrimersMultiplex")
    multiplexes = PCR1Multiplex.objects.using(db_alias).all()
    for mpx in bar(multiplexes):
        old_mpx = OldPrimersMultiplex.objects.using(db_alias).create(name=mpx.name) 
        # This saves selecting the te.
        old_te_ids = [ter.old_te_id for ter in mpx.ters.all()]
        # Batching reduces to selecting the new m2m from the db once.
        old_mpx.primers.add(*old_te_ids)
        transfer_physical_locations(mpx, old_mpx, apps, schema_editor)

@transaction.atomic
def rejoin_panels(apps, schema_editor):
    print()
    print("Reverting panels:")
    db_alias = schema_editor.connection.alias
    Panel = apps.get_model("multiplexes","Panel")
    OldPanel = apps.get_model("linapp","Panel")
    panels = Panel.objects.using(db_alias).all()
    for panel in bar(panels):
        old_panel = OldPanel.objects.using(db_alias).create(name=panel.name) 
        # This saves selecting the te.
        old_te_ids = [te.old_te_id for te in panel.tes.all()]
        # Batching reduces to selecting the new m2m from the db once.
        old_panel.targets.add(*old_te_ids)
