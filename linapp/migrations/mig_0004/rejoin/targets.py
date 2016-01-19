# -*- coding: utf-8 -*-
from frogress import bar

from django.db import IntegrityError, transaction

from linapp.migrations.mig_0004.rejoin.common import unpack_slice, target_types
from linapp.migrations.mig_0004.utils import check_partition, NotCovering, \
    NotMutuallyExclusive, get_partner_ids

def fix_target(target, old_target, apps, schema_editor):
    """
    Perform some common actions for finalizing the transfer:
    Populate the old_target fk; transfer the partners.
    """
    old_target.partner.add(*get_partner_ids(target, apps, schema_editor))
    target.old_target = old_target
    target.save()

def unpack_target(target, apps, schema_editor):
    """
    Create a dict from a target.
    """
    d = unpack_slice(target.slice, apps, schema_editor)
    d["name"] = target.name
    return d

@transaction.atomic
def revert_microsatellites(qs, apps, schema_editor):
    print
    print "Reverting Microsatellites:"
    db_alias = schema_editor.connection.alias
    OldMicrosatellite = apps.get_model("linapp", "Microsatellite")
    type = target_types["Microsatellite",apps,schema_editor]
    for target in bar(qs.select_related("microsatellite")):
        ms = target.microsatellite
        d = unpack_target(target, apps, schema_editor) 
        d.update({k: ms.__dict__[k] for k in [
            "repeat_unit_len",
            "repeat_unit_type",
            "repeat_number",
        ]})
        d["type"] = type
        old_ms = OldMicrosatellite.objects.using(db_alias).create(**d)
        fix_target(target, old_ms, apps, schema_editor)

@transaction.atomic
def revert_known_snps(qs, apps, schema_editor):
    print
    print "Reverting known SNPs:"
    db_alias = schema_editor.connection.alias
    OldSNP = apps.get_model("linapp", "SNP")
    type = target_types["SNP",apps,schema_editor]
    for target in bar(qs.select_related("snp")):
        snp = target.snp
        d = unpack_target(target, apps, schema_editor) 
        d.update({k: snp.__dict__[k] for k in [
            "mutation",
            "modified",
        ]})
        d["type"] = type
        old_snp = OldSNP.objects.using(db_alias).create(**d)
        fix_target(target, old_snp, apps, schema_editor)

def revert_targets_inner(qs, type, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    OldTarget = apps.get_model("linapp", "Target")
    for target in bar(qs):
        d = unpack_target(target, apps, schema_editor)
        d["type"] = type
        old_target = OldTarget.objects.using(db_alias).create(**d)
        fix_target(target, old_target, apps, schema_editor)

@transaction.atomic
def revert_unknown_snps(qs, apps, schema_editor):
    print
    print "Reverting unknown SNPs:"
    type = target_types["SNP",apps,schema_editor]
    revert_targets_inner(qs, type, apps, schema_editor)

@transaction.atomic
def revert_other_targets(qs, apps, schema_editor):
    print
    print "Reverting other targets:"
    # This generates the last type, which we don't use but still want to
    # create.
    target_types["Flank",apps,schema_editor]
    type = target_types["NoSeq",apps,schema_editor]
    revert_targets_inner(qs, type, apps, schema_editor)

def rejoin_targets(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    # We divide into cases. These queries should be a partition of all old
    # targets (mutually exclusive and covering). Otherwise, we fail.
    # FIXME: fix exclude on django Q objects.
    # These funcs fetch a single and multiple target types, respectively.
    Target = apps.get_model("planning", "Target")
    queries = [
        (
            {
                "microsatellite__isnull": False,
                "snp__isnull": True,
            },
            revert_microsatellites,
        ),
        (
            {
                "microsatellite__isnull": True,
                "snp__isnull": False,
                "snp__mutation__isnull": False,
                "snp__modified__isnull": False,
            },
            revert_known_snps,
        ),
        (
            {
                "microsatellite__isnull": True,
                "snp__isnull": False,
                "snp__mutation__isnull": True,
                "snp__modified__isnull": True,
            },
            revert_unknown_snps,
        ),
        (
            {
                "microsatellite__isnull": True,
                "snp__isnull": True,
            },
            revert_other_targets,
        ),
    ]
    targets = Target.objects.using(db_alias).select_related(
        "slice",
        "slice__chromosome",
    )
    try:
        check_partition(targets,[q for q,f in queries])
    except NotCovering:
        raise IntegrityError("We have a wild Target.")
    except NotMutuallyExclusive:
        raise IntegrityError("We have an ambiguous Target.")
    for q,f in queries:
        f(targets.filter(**q), apps, schema_editor)

