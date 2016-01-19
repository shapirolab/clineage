# -*- coding: utf-8 -*-
from frogress import bar

from django.db import transaction

from linapp.migrations.mig_0004.rejoin.common import unpack_slice, target_types

@transaction.atomic
def rejoin_restriction_sites(apps, schema_editor):
    print
    print "Reverting Restriction Sites:"
    db_alias = schema_editor.connection.alias
    # first, migrate the restriction site types.
    RestrictionEnzyme = apps.get_model("planning", "RestrictionEnzyme")
    OldRestrictionSiteType = apps.get_model("linapp", "RestrictionSiteType")
    for renz in RestrictionEnzyme.objects.using(db_alias).all():
        OldRestrictionSiteType.objects.using(db_alias).create(
            **{k: renz.__dict__[k] for k in [
                "id",
                "name",
                "sequence",
                "cut_delta",
                "sticky_bases",
                "sequence_len",
        ]})
    RestrictionSite = apps.get_model("planning", "RestrictionSite")
    OldRestrictionSite = apps.get_model("linapp", "RestrictionSite")
    restriction_sites = RestrictionSite.objects.using(db_alias).select_related(
            "slice",
            "slice__chromosome",
    )
    type = target_types["Plain",apps,schema_editor]
    for restrictionsite in bar(restriction_sites):
        d = unpack_slice(restrictionsite.slice, apps, schema_editor) 
        # Saves need to select_related.
        d["restriction_type_id"] = restrictionsite.enzyme_id
        d["type"] = type
        old_restrictionsite = OldRestrictionSite.objects.using(db_alias) \
            .create(**d)
