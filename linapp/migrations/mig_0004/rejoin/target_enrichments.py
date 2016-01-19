# -*- coding: utf-8 -*-
from frogress import bar

from django.db import transaction, IntegrityError

from linapp.migrations.mig_0004.utils import transfer_physical_locations, \
    getdna, get_partner_ids

def create_target_enrichment(ter, type, apps, schema_editor):
    """
    Create an old target enrichment, incorporating a new ter and te.
    """
    db_alias = schema_editor.connection.alias
    OldTargetEnrichment = apps.get_model("linapp", "TargetEnrichment")
    old_te = OldTargetEnrichment.objects.using(db_alias).create(
        type=type,
        amplicon=getdna(
            ter.te.chromosome,
            ter.te.left.slice.start_pos,
            ter.te.right.slice.end_pos,
        ),
        chromosome_id=ter.te.chromosome_id,
        left_id=ter.left_primer.old_primer_id,
        right_id=ter.right_primer.old_primer_id,
        **{k: ter.__dict__[k] for k in [
            "passed_validation",
            "validation_failure_id",
            "validation_date",
            "comment",
        ]}
    )
    target_ids = [target.old_target_id for target in \
        ter.te.targets.using(db_alias).all()]
    old_te.targets.add(*target_ids)
    old_te.partner.add(*get_partner_ids(ter.te, apps, schema_editor))
    ter.old_te = old_te
    ter.save()
    # To save select on old_te.
    if ter.te.old_te_id is None:
        ter.te.old_te = old_te
        ter.te.save()
    transfer_physical_locations(ter, old_te, apps, schema_editor)
    return old_te

@transaction.atomic
def revert_ters(qs, type_name, apps, schema_editor):
    """
    Revert all TERs in a given qs, given their type name.
    """
    db_alias = schema_editor.connection.alias
    OldTargetEnrichmentType = apps.get_model("linapp", "TargetEnrichmentType")
    type, c = OldTargetEnrichmentType.objects.using(db_alias).get_or_create(
        name=type_name,
    )
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
        create_target_enrichment(ter, type, apps, schema_editor)

def revert_pcr1_primer_pair_ter(apps, schema_editor):
    print
    print "Reverting PCR1 Primer Pair TERs:"
    db_alias = schema_editor.connection.alias
    PCR1PrimerPairTER = apps.get_model("reagents", "PCR1PrimerPairTER")
    # Here we divide into cases, but there are none.
    qs = PCR1PrimerPairTER.objects.using(db_alias).all()
    revert_ters(qs, "PCR_with_tails", apps, schema_editor)

def revert_pcr1_with_company_tag_primer_pair_ter(apps, schema_editor):
    print
    print "Reverting PCR1 With Company Tag Primer Pair TERs:"
    db_alias = schema_editor.connection.alias
    PCR1WithCompanyTagPrimerPairTER = apps.get_model("reagents", "PCR1WithCompanyTagPrimerPairTER")
    # Here we divide into cases, but there are none.
    qs = PCR1WithCompanyTagPrimerPairTER.objects.using(db_alias).all()
    revert_ters(qs, "PCR_with_tails", apps, schema_editor)

def revert_pcr1_primer_pair_ter_deprecated(apps, schema_editor):
    print
    print "Reverting Deprecated PCR1 Primer Pair TERs:"
    db_alias = schema_editor.connection.alias
    PCR1PrimerPairTERDeprecated = apps.get_model("reagents", "PCR1PrimerPairTERDeprecated")
    # Here we divide into cases, but there are none.
    qs = PCR1PrimerPairTERDeprecated.objects.using(db_alias).all()
    revert_ters(qs, "depricated", apps, schema_editor)  # sic.

def revert_no_tail_ter(apps, schema_editor):
    print
    print "Reverting No-Tails Primer Pair TERs:"
    db_alias = schema_editor.connection.alias
    TargetedNoTailPrimerPairTER = apps.get_model("reagents", "TargetedNoTailPrimerPairTER")
    # Here we divide into cases, but there are none.
    qs = TargetedNoTailPrimerPairTER.objects.using(db_alias).all()
    revert_ters(qs, "PCR", apps, schema_editor)

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
