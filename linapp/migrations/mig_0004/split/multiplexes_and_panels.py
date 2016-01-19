# -*- coding: utf-8 -*-
from frogress import bar

from django.db import IntegrityError, transaction

from linapp.migrations.mig_0004.utils import transfer_physical_locations

@transaction.atomic
def split_multiplexes(apps, schema_editor):
    print
    print "Converting PCR1 multiplexes:"
    db_alias = schema_editor.connection.alias
    OldPrimersMultiplex = apps.get_model("linapp","PrimersMultiplex")
    PCR1Multiplex = apps.get_model("multiplexes","PCR1Multiplex")
    old_multiplexes = OldPrimersMultiplex.objects.using(db_alias).all()
    for old_mpx in bar(old_multiplexes):
        mpx = PCR1Multiplex.objects.using(db_alias).create(name=old_mpx.name) 
        q = {
            "new_ter_model__in": [
                "PCR1PrimerPairTER",
                "PCR1WithCompanyTagPrimerPairTER",
                "PCR1PrimerPairTERDeprecated",
            ],
        }
        bads = old_mpx.primers.exclude(**q)
        if bads:
            raise IntegrityError("PCR1 Multiplex has bad TERs, ex. {}".format(
                bads[0].id))
        # This saves selecting the ter.
        ter_ids = [old_te.new_ter_id for old_te in old_mpx.primers.filter(**q)]
        # Batching reduces to selecting the new m2m from the db once.
        mpx.ters.add(*ter_ids)
        transfer_physical_locations(old_mpx, mpx, apps, schema_editor)

@transaction.atomic
def split_panels(apps, schema_editor):
    print
    print "Converting panels:"
    db_alias = schema_editor.connection.alias
    OldPanel = apps.get_model("linapp","Panel")
    Panel = apps.get_model("multiplexes","Panel")
    old_panels = OldPanel.objects.using(db_alias).all()
    for old_panel in bar(old_panels):
        panel = Panel.objects.using(db_alias).create(name=old_panel.name) 
        # This saves selecting the te.
        te_ids = [old_te.new_te_id for old_te in old_panel.targets.all()]
        # Batching reduces to selecting the new m2m from the db once.
        panel.tes.add(*te_ids)
