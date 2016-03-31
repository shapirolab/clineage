# -*- coding: utf-8 -*-
from frogress import bar

from django.db import transaction, IntegrityError

from linapp.migrations.mig_0004.utils import transfer_physical_locations, \
    getdna, get_partner_ids, rc
from linapp.migrations.mig_0004.rejoin.common import unpack_slice, \
    get_or_create_sequence, target_types, FuncCacheDict
# NOTE: from various other places!
from linapp.migrations.mig_0004.split.primers import LEFT_TAIL, RIGHT_TAIL
from linapp.migrations.mig_0004.rejoin.primers import primer_tails

def _get_or_create_target_enrichment_type(type, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    OldTargetEnrichmentType = apps.get_model("linapp","TargetEnrichmentType")
    old_type, c = OldTargetEnrichmentType.objects.using(db_alias).get_or_create(name=type)
    return old_type

te_types = FuncCacheDict(_get_or_create_target_enrichment_type)

def create_target_enrichment_from_ter(ter, type_name, apps, schema_editor):
    """
    Create an old target enrichment, incorporating a new ter and te.
    """
    db_alias = schema_editor.connection.alias
    OldTargetEnrichment = apps.get_model("linapp", "TargetEnrichment")
    d = {k: ter.__dict__[k] for k in [
            "passed_validation",
            "validation_failure_id",
            "validation_date",
            "comment",
    ]}
    if ter.old_adam_te_pk is not None:
        d["id"] = ter.old_adam_te_pk
    old_te = OldTargetEnrichment.objects.using(db_alias).create(
        type=te_types[type_name,apps,schema_editor],
        amplicon=getdna(
            ter.te.chromosome,
            ter.te.left.slice.start_pos,
            ter.te.right.slice.end_pos,
        ),
        chromosome_id=ter.te.chromosome_id,
        left_id=ter.left_primer.old_primer_id,
        right_id=ter.right_primer.old_primer_id,
        **d
    )
    target_ids = [target.old_target_id for target in \
        ter.te.targets.using(db_alias).all()]
    old_te.targets.add(*target_ids)
    old_te.partner.add(*get_partner_ids(ter.te, apps, schema_editor))
    ter.old_te = old_te
    ter.save()
    # To save select on old_te.
    if ter.te.old_te_id is None or ter.old_adam_te_pk is not None:
        ter.te.old_te = old_te
        ter.te.save()
    transfer_physical_locations(ter, old_te, apps, schema_editor)
    return old_te

def create_target_enrichment_from_te(te, apps, schema_editor):
    """
    Create an old target enrichment out of a te with no ter.
    Requires creating primers from the UGSs, has no physical location.
    """
    db_alias = schema_editor.connection.alias
    # We guess which type of te this was.
    if te.planning_version == 0:
        # This means it was planned in the no-tails era.
        fwd_tail = ""
        rev_tail = ""
        type = te_types["PCR",apps,schema_editor]
    else:
        # This means it was planned in the PCR1 (with-tails) era.
        # NOTE: this might apply to new types of tes as well.
        fwd_tail = LEFT_TAIL
        rev_tail = RIGHT_TAIL
        type = te_types["PCR_with_tails",app,schema_editor]
    # First we generate the primers from the UGSs.
    # No physical locations!
    OldPrimer = apps.get_model("linapp","Primer")
    d = unpack_slice(te.left.slice, apps, schema_editor)
    head = d["referencevalue"].sequence
    d["name"] = "Mig_auto_fwd_{}".format(te.id)
    d["type"] = target_types["Plain", apps, schema_editor]
    d["strand"] = "+"
    d["sequence"] = get_or_create_sequence(fwd_tail+head, apps, schema_editor)
    d["tail"] = primer_tails[fwd_tail,apps,schema_editor]
    old_left = OldPrimer.objects.using(db_alias).create(**d)
    d = unpack_slice(te.right.slice, apps, schema_editor)
    head = rc(d["referencevalue"].sequence)
    d["name"] = "Mig_auto_rev_{}".format(te.id)
    d["type"] = target_types["Plain", apps, schema_editor]
    d["strand"] = "-"
    d["sequence"] = get_or_create_sequence(rev_tail+head, apps, schema_editor)
    d["tail"] = primer_tails[rev_tail,apps,schema_editor]
    old_right = OldPrimer.objects.using(db_alias).create(**d)
    # Now we can generate the te.
    OldTargetEnrichment = apps.get_model("linapp", "TargetEnrichment")
    old_te = OldTargetEnrichment.objects.using(db_alias).create(
        type=type,
        amplicon=getdna(
            te.chromosome,
            te.left.slice.start_pos,
            te.right.slice.end_pos,
        ),
        chromosome_id=te.chromosome_id,
        left=old_left,
        right=old_right,
    )
    target_ids = [target.old_target_id for target in \
        te.targets.using(db_alias).all()]
    old_te.targets.add(*target_ids)
    old_te.partner.add(*get_partner_ids(te, apps, schema_editor))
    # No physical locations!
    te.old_te = old_te
    te.save()
    return old_te

@transaction.atomic
def revert_ters(qs, type_name, apps, schema_editor):
    """
    Revert all TERs in a given qs, given their type name.
    """
    db_alias = schema_editor.connection.alias
    ters = qs.select_related(
        "te",
        "te__chromosome",
        "te__left",
        "te__left__slice",
        "te__right",
        "te__right__slice",
        "left_primer",
        "right_primer",
    )
    for ter in bar(ters):
        create_target_enrichment_from_ter(ter, type_name, apps, schema_editor)

def revert_pcr1_primer_pair_ter(apps, schema_editor):
    print()
    print("Reverting PCR1 Primer Pair TERs:")
    db_alias = schema_editor.connection.alias
    PCR1PrimerPairTER = apps.get_model("reagents", "PCR1PrimerPairTER")
    # Here we divide into cases, but there are none.
    qs = PCR1PrimerPairTER.objects.using(db_alias).all()
    revert_ters(qs, "PCR_with_tails", apps, schema_editor)

def revert_pcr1_with_company_tag_primer_pair_ter(apps, schema_editor):
    print()
    print("Reverting PCR1 With Company Tag Primer Pair TERs:")
    db_alias = schema_editor.connection.alias
    PCR1WithCompanyTagPrimerPairTER = apps.get_model("reagents", "PCR1WithCompanyTagPrimerPairTER")
    # Here we divide into cases, but there are none.
    qs = PCR1WithCompanyTagPrimerPairTER.objects.using(db_alias).all()
    revert_ters(qs, "PCR_with_tails", apps, schema_editor)

def revert_pcr1_primer_pair_ter_deprecated(apps, schema_editor):
    print()
    print("Reverting Deprecated PCR1 Primer Pair TERs:")
    db_alias = schema_editor.connection.alias
    PCR1PrimerPairTERDeprecated = apps.get_model("reagents", "PCR1PrimerPairTERDeprecated")
    # Here we divide into cases, but there are none.
    qs = PCR1PrimerPairTERDeprecated.objects.using(db_alias).all()
    revert_ters(qs, "depricated", apps, schema_editor)  # sic.

def revert_no_tail_ter(apps, schema_editor):
    print()
    print("Reverting No-Tails Primer Pair TERs:")
    db_alias = schema_editor.connection.alias
    TargetedNoTailPrimerPairTER = apps.get_model("reagents", "TargetedNoTailPrimerPairTER")
    # Here we divide into cases, but there are none.
    qs = TargetedNoTailPrimerPairTER.objects.using(db_alias).all()
    revert_ters(qs, "PCR", apps, schema_editor)

@transaction.atomic
def revert_tes_without_ter(apps, schema_editor):
    """
    Revert all TEs which don't have a TER (so they weren't reverted before).
    """
    print()
    print("Reverting TEs without TERs:")
    db_alias = schema_editor.connection.alias
    TargetEnrichment = apps.get_model("planning", "TargetEnrichment")
    # We only take those that weren't assigned an old_te, meaning they
    # didn't have a TER and therefore weren't reverted.
    tes = TargetEnrichment.objects.using(db_alias).filter(old_te__isnull=True).select_related(
        "chromosome",
        "left",
        "left__slice",
        "left__slice__chromosome",
        "right",
        "right__slice",
        "right__slice__chromosome",
    )
    for te in bar(tes):
        create_target_enrichment_from_te(te, apps, schema_editor)

def rejoin_target_enrichments(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    # first, migrate the restriction site types.
    TargetEnrichmentFailureType = apps.get_model("reagents", "TargetEnrichmentFailureType")
    OldTargetEnrichmentFailureType = apps.get_model("linapp", "TargetEnrichmentFailureType")
    with transaction.atomic():
        for teft in TargetEnrichmentFailureType.objects.using(db_alias).all():
            OldTargetEnrichmentFailureType.objects.using(db_alias).create(
                **{k: teft.__dict__[k] for k in [
                    "id",
                    "name",
                    "description",
            ]})
    # Check for stray TERs.
    PCR1PrimerPairTERBase = apps.get_model("reagents", "PCR1PrimerPairTERBase")
    if PCR1PrimerPairTERBase.objects.using(db_alias).filter(
        pcr1primerpairter__isnull=True,
        pcr1withcompanytagprimerpairter__isnull=True,
        pcr1primerpairterdeprecated__isnull=True,
    ):
        raise IntegrityError("We have a stray PCR1PrimerPairTERBase.")
    revert_pcr1_primer_pair_ter(apps, schema_editor)
    revert_pcr1_with_company_tag_primer_pair_ter(apps, schema_editor)
    revert_pcr1_primer_pair_ter_deprecated(apps, schema_editor)
    revert_no_tail_ter(apps, schema_editor)
    revert_tes_without_ter(apps, schema_editor)
