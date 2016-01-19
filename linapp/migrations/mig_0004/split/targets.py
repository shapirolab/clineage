# -*- coding: utf-8 -*-
from frogress import bar

from django.db import IntegrityError, transaction
from django.db.models import F

from linapp.migrations.mig_0004.split.common import get_or_create_slice
from linapp.migrations.mig_0004.split.primers import split_primers
from linapp.migrations.mig_0004.utils import check_partition, NotCovering, \
    NotMutuallyExclusive, get_partner_ids

def fix_target(old_target, target, apps, schema_editor):
    """
    Perform some common actions for finalizing the transfer:
    Populate the new_target fk; transfer the partners.
    """
    target.partner.add(*get_partner_ids(old_target, apps, schema_editor))
    old_target.new_target = target
    old_target.save()

def prepare_target_dict(old_target, apps, schema_editor):
    """
    Prepare the common fields for a new target.
    """
    slice = get_or_create_slice(old_target, apps, schema_editor)
    return {
        "name": old_target.name,
        "slice": slice,
    }

@transaction.atomic
def convert_microsatellites(qs, apps, schema_editor):
    print
    print "Converting Microsatellites:"
    db_alias = schema_editor.connection.alias
    Microsatellite = apps.get_model("planning", "Microsatellite")
    for old_target in bar(qs.select_related("microsatellite")):
        old_ms = old_target.microsatellite
        d = prepare_target_dict(old_target, apps, schema_editor) 
        d.update({k: old_ms.__dict__[k] for k in [
            "repeat_unit_len",
            "repeat_unit_type",
            "repeat_number",
        ]})
        ms = Microsatellite.objects.using(db_alias).create(**d)
        fix_target(old_target, ms, apps, schema_editor)

@transaction.atomic
def convert_known_snps(qs, apps, schema_editor):
    print
    print "Converting known SNPs:"
    db_alias = schema_editor.connection.alias
    SNP = apps.get_model("planning", "SNP")
    for old_target in bar(qs.select_related("snp")):
        old_snp = old_target.snp
        d = prepare_target_dict(old_target, apps, schema_editor) 
        d.update({k: old_snp.__dict__[k] for k in [
            "mutation",
            "modified",
        ]})
        snp = SNP.objects.using(db_alias).create(**d)
        fix_target(old_target, snp, apps, schema_editor)

@transaction.atomic
def convert_unknown_snps(qs, apps, schema_editor):
    print
    print "Converting unknown SNPs:"
    db_alias = schema_editor.connection.alias
    SNP = apps.get_model("planning", "SNP")
    for old_target in bar(qs):
        d = prepare_target_dict(old_target, apps, schema_editor) 
        snp = SNP.objects.using(db_alias).create(**d)
        fix_target(old_target, snp, apps, schema_editor)

@transaction.atomic
def convert_restrictionsites(qs, apps, schema_editor):
    print
    print "Converting Restriction Sites:"
    # first, migrate the restriction site types.
    db_alias = schema_editor.connection.alias
    OldRestrictionSiteType = apps.get_model("linapp", "RestrictionSiteType")
    RestrictionEnzyme = apps.get_model("planning", "RestrictionEnzyme")
    for orst in OldRestrictionSiteType.objects.using(db_alias).all():
        RestrictionEnzyme.objects.using(db_alias).create(
            **{k: orst.__dict__[k] for k in [
                "id",
                "name",
                "sequence",
                "cut_delta",
                "sticky_bases",
                "sequence_len",
        ]})
    RestrictionSite = apps.get_model("planning", "RestrictionSite")
    for old_target in bar(qs.select_related("restrictionsite")):
        old_restrictionsite = old_target.restrictionsite
        slice = get_or_create_slice(old_target, apps, schema_editor) 
        restrictionsite = RestrictionSite.objects.using(db_alias).create(
            slice=slice,
            # No select_related required, we only need the id.
            enzyme_id=old_restrictionsite.restriction_type_id,
        )

def split_targets(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    # We divide into cases. These queries should be a partition of all old
    # targets (mutually exclusive and covering). Otherwise, we fail.
    # FIXME: fix exclude on django Q objects.
    # These funcs fetch a single and multiple target types, respectively.
    OldTarget = apps.get_model("linapp", "Target")
    OldTargetType = apps.get_model("linapp", "TargetType")
    def get_tt(x):
        return OldTargetType.objects.using(db_alias).get(name=x)
    def get_tts(xs):
        return [OldTargetType.objects.using(db_alias).get(name=x) for x in xs]
    queries = [
        (
            {
                "microsatellite__isnull": False,
                "snp__isnull": True,
                "restrictionsite__isnull": True,
                "primer__isnull": True,
                "type": get_tt("Microsatellite"),
            },
            convert_microsatellites,
        ),
        (
            {
                "microsatellite__isnull": True,
                "snp__isnull": False,
                "restrictionsite__isnull": True,
                "primer__isnull": True,
                "type": get_tt("SNP"),
            },
            convert_known_snps,
        ),
        (
            {
                "microsatellite__isnull": True,
                "snp__isnull": True,
                "restrictionsite__isnull": True,
                "primer__isnull": True,
                "type__in": get_tts(["SNP","NoSeq"]),
                "start_pos": F("end_pos"),
            },
            convert_unknown_snps,
        ),
        (
            {
                "microsatellite__isnull": True,
                "snp__isnull": True,
                "restrictionsite__isnull": False,
                "primer__isnull": True,
                "type": get_tt("Plain"),
                },
            convert_restrictionsites,
        ),
        (
            {
                "microsatellite__isnull": True,
                "snp__isnull": True,
                "restrictionsite__isnull": True,
                "primer__isnull": False,
                "type__in": get_tts(["Plain","Flank"]),
            },
            split_primers,
        ),
    ]
    old_targets = OldTarget.objects.using(db_alias)
    try:
        check_partition(old_targets,[q for q,f in queries])
    except NotCovering as e:
        raise IntegrityError("We have a wild Target, example: {}".format(
            e.args[0][0].id))
    except NotMutuallyExclusive as e:
        raise IntegrityError("We have an ambiguous Target - queries #{} and " \
            "#{} both match. Example: {}".format(e.args[0],e.args[1],
                e.args[2][0].id))
    for q,f in queries:
        f(old_targets.filter(**q), apps, schema_editor)
