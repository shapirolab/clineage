# -*- coding: utf-8 -*-
from frogress import bar

from django.db import IntegrityError, transaction
from django.db.models import Count

from linapp.migrations.mig_0004.rejoin.common import unpack_slice, \
    get_or_create_sequence, target_types, FuncCacheDict
from linapp.migrations.mig_0004.utils import transfer_physical_locations, rc

SELECT_RELATED_NOTAIL = [
    "ugs",
    "ugs__slice",
    "ugs__slice__chromosome",
]

SELECT_RELATED_TAIL = SELECT_RELATED_NOTAIL + [
    "irac",
    "irac__ira",
]
    
def unpack_common_primer(primer, apps, schema_editor):
    """
    Create a dict from a primer.
    """
    d = unpack_slice(primer.ugs.slice, apps, schema_editor)
    d["name"] = primer.name
    d["type"] = target_types["Plain", apps, schema_editor]
    return d

def _get_or_create_primer_tail(tail, apps, schema_editor):
    if tail == "":
        return None
    db_alias = schema_editor.connection.alias
    OldPrimerTail = apps.get_model("linapp","PrimerTail")
    old_tail, c = OldPrimerTail.objects.using(db_alias).get_or_create(tail=tail)
    return old_tail

primer_tails = FuncCacheDict(_get_or_create_primer_tail)

def create_fwd_primer(primer, apps, schema_editor, tail="", extra=""):
    """
    Create an old forward primer, given the actual sequence of the tail, and
    taking the head from the reference value of the target.  If extra is given,
    it is added to the sequence but not to the tail.
    Also transfers the locations.
    """
    db_alias = schema_editor.connection.alias
    OldPrimer = apps.get_model("linapp","Primer")
    d = unpack_common_primer(primer, apps, schema_editor)
    head = d["referencevalue"].sequence
    d["strand"] = "+"
    d["sequence"] = get_or_create_sequence(tail+extra+head, apps, schema_editor)
    d["tail"] = primer_tails[tail,apps,schema_editor]
    old_primer = OldPrimer.objects.using(db_alias).create(**d)
    primer.old_primer = old_primer
    primer.save()
    # To save select on old_primer.
    if primer.ugs.old_primer_id is None:
        primer.ugs.old_primer = old_primer
        primer.ugs.save()
    transfer_physical_locations(primer, old_primer, apps, schema_editor)

def create_rev_primer(primer, apps, schema_editor, tail="", extra=""):
    """
    Create an old reverse primer, given the actual sequence of the tail, and
    taking the head from the reference value of the target.  If extra is given,
    it is added to the sequence but not to the tail.
    Also transfers the locations.
    """
    db_alias = schema_editor.connection.alias
    OldPrimer = apps.get_model("linapp","Primer")
    d = unpack_common_primer(primer, apps, schema_editor)
    head = rc(d["referencevalue"].sequence)
    d["strand"] = "-"
    d["sequence"] = get_or_create_sequence(tail+extra+head, apps, schema_editor)
    d["tail"] = primer_tails[tail,apps,schema_editor]
    old_primer = OldPrimer.objects.using(db_alias).create(**d)
    primer.old_primer = old_primer
    primer.save()
    # To save select on old_primer.
    if primer.ugs.old_primer_id is None:
        primer.ugs.old_primer = old_primer
        primer.ugs.save()
    transfer_physical_locations(primer, old_primer, apps, schema_editor)

def get_fwd_irac_tail(primer):
    """
    Get the pcr1 tail from a the irac of a forward primer.
    """
    return primer.irac.ira._sequence[primer.irac.overlap_start:]

def get_rev_irac_tail(primer):
    """
    Get the pcr1 tail from a the irac of a reverse primer.
    """
    return rc(primer.irac.ira._sequence)[primer.irac.overlap_start:]

@transaction.atomic
def revert_pcr1_primers(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    PCR1PlusPrimer = apps.get_model("synthesis","PCR1PlusPrimer")
    PCR1MinusPrimer = apps.get_model("synthesis","PCR1MinusPrimer")
    plus_qs = PCR1PlusPrimer.objects.using(db_alias).select_related(
        *SELECT_RELATED_TAIL)
    minus_qs = PCR1MinusPrimer.objects.using(db_alias).select_related(
        *SELECT_RELATED_TAIL)
    print
    print "Reverting forward PCR1 primers:"
    for primer in bar(plus_qs):
        tail = get_fwd_irac_tail(primer)
        create_fwd_primer(primer, apps, schema_editor, tail)
    print
    print "Reverting reverse PCR1 primers:"
    for primer in bar(minus_qs):
        tail = get_rev_irac_tail(primer)
        create_rev_primer(primer, apps, schema_editor, tail)

@transaction.atomic
def revert_pcr1_with_company_tag_primers(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    PCR1WithCompanyTagPlusPrimer = apps.get_model("synthesis","PCR1WithCompanyTagPlusPrimer")
    PCR1WithCompanyTagMinusPrimer = apps.get_model("synthesis","PCR1WithCompanyTagMinusPrimer")
    plus_qs = PCR1WithCompanyTagPlusPrimer.objects.using(db_alias) \
        .select_related(*SELECT_RELATED_TAIL)
    minus_qs = PCR1WithCompanyTagMinusPrimer.objects.using(db_alias) \
        .select_related(*SELECT_RELATED_TAIL)
    print
    print "Reverting forward PCR1 primers with company tags:"
    for primer in bar(plus_qs):
        tail = get_fwd_irac_tail(primer)
        create_fwd_primer(primer, apps, schema_editor, tail, extra=primer.tag)
    print
    print "Reverting reverse PCR1 primers with company tags:"
    for primer in bar(minus_qs):
        tail = get_rev_irac_tail(primer)
        create_rev_primer(primer, apps, schema_editor, tail, extra=primer.tag)

@transaction.atomic
def revert_no_tail_primers(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    TargetedNoTailPlusPrimer = apps.get_model("synthesis","TargetedNoTailPlusPrimer")
    TargetedNoTailMinusPrimer = apps.get_model("synthesis","TargetedNoTailMinusPrimer")
    plus_qs = TargetedNoTailPlusPrimer.objects.using(db_alias) \
        .select_related(*SELECT_RELATED_NOTAIL)
    minus_qs = TargetedNoTailMinusPrimer.objects.using(db_alias) \
        .select_related(*SELECT_RELATED_NOTAIL)
    print
    print "Reverting forward PCR1 primers:"
    for primer in bar(plus_qs):
        create_fwd_primer(primer, apps, schema_editor)
    print
    print "Reverting reverse PCR1 primers:"
    for primer in bar(minus_qs):
        create_rev_primer(primer, apps, schema_editor)

def validate_ugss_without_primers(apps, schema_editor):
    """
    Validate the UGSs which don't have an old_primer, asserting they have a TE.
    This is required so that the TE will revert them to Primers later.
    """
    db_alias = schema_editor.connection.alias
    UGSPlus = apps.get_model("planning","UGSPlus")
    UGSMinus = apps.get_model("planning","UGSMinus")
    plus_qs = UGSPlus.objects.using(db_alias) \
        .filter(old_primer__isnull=True) \
        .select_related("slice","slice__chromosome") \
        .annotate(te_c=Count("targetenrichment")).filter(te_c=0)
    minus_qs = UGSMinus.objects.using(db_alias) \
        .filter(old_primer__isnull=True) \
        .select_related("slice","slice__chromosome") \
        .annotate(te_c=Count("targetenrichment")).filter(te_c=0)
    print
    print "Validating forward UGSs without primers."
    if plus_qs:
        raise IntegrityError("Bad UGS - no Primer, no TE. ex. {}".format(
            plus_qs[0].id))
    print
    print "Validating reverse UGSs without primers."
    if minus_qs:
        raise IntegrityError("Bad UGS - no Primer, no TE. ex. {}".format(
            minus_qs[0].id))

def rejoin_primers(apps, schema_editor):
    revert_pcr1_primers(apps, schema_editor)
    revert_pcr1_with_company_tag_primers(apps, schema_editor)
    revert_no_tail_primers(apps, schema_editor)
    validate_ugss_without_primers(apps, schema_editor)
