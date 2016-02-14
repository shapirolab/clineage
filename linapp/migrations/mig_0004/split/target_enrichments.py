# -*- coding: utf-8 -*-
from frogress import bar

from django.db import IntegrityError, transaction
from django.db.models import F

from linapp.migrations.mig_0004.utils import get_physical_locations, \
    alter_physical_locations, check_partition, NotCovering, \
    NotMutuallyExclusive, get_partner_ids

def create_te_and_ter(old_te, ter_model, planning_version, apps, schema_editor, force_loc=False):
    db_alias = schema_editor.connection.alias
    TargetEnrichment = apps.get_model("planning", "TargetEnrichment")
    old_left = old_te.left
    old_right = old_te.right
    te, c = TargetEnrichment.objects.using(db_alias).get_or_create(
        chromosome_id=old_te.chromosome_id,
        left_id=old_left.new_ugs_id,
        right_id=old_right.new_ugs_id,
        planning_version=planning_version,
    )
    target_ids = [old_target.new_target_id for old_target in \
        old_te.targets.using(db_alias).all()]
    te.targets.add(*target_ids)
    te.partner.add(*get_partner_ids(old_te, apps, schema_editor))
    old_te.new_te = te
    ret = te
    locations = get_physical_locations(old_te, apps, schema_editor)
    if locations:
        if ter_model is None:
            raise IntegrityError("We have a TE with a physical_location, but"
                "whose primers don't have. ex. {}".format(old_te.id))
        TargetEnrichmentReagent = apps.get_model("reagents", ter_model)
        ter = TargetEnrichmentReagent.objects.using(db_alias).create(
            te=te,
            left_primer_id=old_left.new_primer_id,
            right_primer_id=old_right.new_primer_id,
            old_adam_te_pk=old_te.id,  # Temporary field to assist with the migration of the adamiya.
            **{k: old_te.__dict__[k] for k in [
                "passed_validation",
                # Assumes the TargetEnrichmentFailureType-s have already been
                # migrated, preserving their id.
                "validation_failure_id",
                "validation_date",
                "comment",
            ]}
        )
        ret = ter
        old_te.new_ter_id = ter.id
        old_te.new_ter_model = ter_model
        alter_physical_locations(ter, locations, apps, schema_editor)
    elif force_loc:
        raise IntegrityError("We have a bad TE without a physical_location, "
            "ex. {}".format(old_te.id))
    old_te.save()
    return ret

@transaction.atomic
def convert_pcr1_target_enrichments(qs, apps, schema_editor):
    print
    print "Converting normal PCR1 target enrichments:"
    for old_te in bar(qs):
        create_te_and_ter(old_te, "PCR1PrimerPairTER", 1, apps, schema_editor)

@transaction.atomic
def convert_pcr1_with_tag_target_enrichments(qs, apps, schema_editor):
    print
    print "Converting PCR1 target enrichments with company tag:"
    for old_te in bar(qs):
        create_te_and_ter(old_te, "PCR1WithCompanyTagPrimerPairTER", 1, apps, schema_editor)

@transaction.atomic
def convert_deprecated_pcr1_target_enrichments(qs, apps, schema_editor):
    print
    print "Converting deprecated PCR1 target enrichments:"
    for old_te in bar(qs):
        create_te_and_ter(old_te, "PCR1PrimerPairTERDeprecated", 1, apps, schema_editor, force_loc=True)

@transaction.atomic
def convert_no_tail_target_enrichments(qs, apps, schema_editor):
    print
    print "Converting no-tail target enrichments:"
    for old_te in bar(qs):
        create_te_and_ter(old_te, "TargetedNoTailPrimerPairTER", 0, apps, schema_editor)

@transaction.atomic
def convert_theoretical_pcr1_target_enrichments(qs, apps, schema_editor):
    print
    print "Converting theoretical PCR1 target enrichments:"
    for old_te in bar(qs):
        create_te_and_ter(old_te, None, 1, apps, schema_editor)

@transaction.atomic
def convert_theoretical_no_tail_target_enrichments(qs, apps, schema_editor):
    print
    print "Converting theoretical no-tail target enrichments:"
    for old_te in bar(qs):
        create_te_and_ter(old_te, None, 0, apps, schema_editor)

def split_target_enrichments(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    # first, migrate the target enrichment failure types.
    OldTargetEnrichmentFailureType  = apps.get_model("linapp", "TargetEnrichmentFailureType")
    TargetEnrichmentFailureType = apps.get_model("reagents", "TargetEnrichmentFailureType")
    with transaction.atomic():
        for oteft in OldTargetEnrichmentFailureType.objects.using(db_alias).all():
            TargetEnrichmentFailureType.objects.using(db_alias).create(
                **{k: oteft.__dict__[k] for k in [
                    "id",
                    "name",
                    "description",
            ]})
    # Now we divide into cases. These queries should be a partition of all old
    # targeted_enrichments (mutually exclusive and covering). Otherwise, we
    # fail.
    # FIXME: fix exclude on django Q objects.
    # This func fetch target enrichment types.
    OldTargetEnrichment = apps.get_model("linapp", "TargetEnrichment")
    OldTargetEnrichmentType = apps.get_model("linapp", "TargetEnrichmentType")
    def get_tet(x):
        return OldTargetEnrichmentType.objects.using(db_alias).get(name=x)
    queries = [
        (
            {
                "left__new_primer_model": "PCR1PlusPrimer",
                "right__new_primer_model": "PCR1MinusPrimer",
                "left__chromosome": F("chromosome"),
                "right__chromosome": F("chromosome"),
                "type": get_tet("PCR_with_tails"),
            },
            convert_pcr1_target_enrichments,
        ),
        (
            {
                "left__new_primer_model": "PCR1WithCompanyTagPlusPrimer",
                "right__new_primer_model": "PCR1WithCompanyTagMinusPrimer",
                "left__chromosome": F("chromosome"),
                "right__chromosome": F("chromosome"),
                "type": get_tet("PCR_with_tails"),
            },
            convert_pcr1_with_tag_target_enrichments,
        ),
        (
            {
                "left__new_primer_model": "PCR1PlusPrimer",
                "right__new_primer_model": "PCR1MinusPrimer",
                "left__chromosome": F("chromosome"),
                "right__chromosome": F("chromosome"),
                "type": get_tet("depricated"),  # sic.
            },
            convert_deprecated_pcr1_target_enrichments,
        ),
        (
            {
                "left__new_primer_model": "TargetedNoTailPlusPrimer",
                "right__new_primer_model": "TargetedNoTailMinusPrimer",
                "left__chromosome": F("chromosome"),
                "right__chromosome": F("chromosome"),
                "type": get_tet("PCR"),
            },
            convert_no_tail_target_enrichments,
        ),
        (
            {
                "left__new_primer_id__isnull": True,
                "right__new_primer_id__isnull": True,
                "left__chromosome": F("chromosome"),
                "right__chromosome": F("chromosome"),
                "type": get_tet("PCR_with_tails"),
            },
            convert_theoretical_pcr1_target_enrichments,
        ),
        (
            {
                "left__new_primer_id__isnull": True,
                "right__new_primer_id__isnull": True,
                "left__chromosome": F("chromosome"),
                "right__chromosome": F("chromosome"),
                "type": get_tet("PCR"),
            },
            convert_theoretical_no_tail_target_enrichments,
        ),
    ]
    old_tes = OldTargetEnrichment.objects.using(db_alias).select_related(
        "left",
        "right",
    )
    try:
        check_partition(old_tes,[q for q,f in queries])
    except NotCovering as e:
        raise IntegrityError("We have a TargetEnrichment with unmatching "
            "primers, example: {}".format(e.args[0][0].id))
    except NotMutuallyExclusive as e:
        raise RuntimeError("WTF? Queries #{} and #{} have things in common, "
            "ex. {}".format(e.args[0],e.args[1],e.args[2][0].id))
    for q,f in queries:
        f(old_tes.filter(**q), apps, schema_editor)

