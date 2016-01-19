# -*- coding: utf-8 -*-
from frogress import bar
import re

from django.db import IntegrityError, transaction
from django.db.models import Count

from linapp.migrations.mig_0004.split.common import get_or_create_slice
from linapp.migrations.mig_0004.utils import transfer_physical_locations, \
    check_partition, NotCovering, NotMutuallyExclusive, getdna, rc

LEFT_TAIL = "CTACACGACGCTCTTCCGATCT" # FIXME!
RIGHT_TAIL = "CAGACGTGTGCTCTTCCGATCT" # FIXME!
IRA1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" # FIXME!
IRA2 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # FIXME!
IRAC1_OVERLAP_START = 11 # FIXME
IRAC1_OVERLAP_END = 33 # FIXME
IRAC2_OVERLAP_START = 12 # FIXME
IRAC2_OVERLAP_END = 34 # FIXME

def check_left_primer(old_primer, tail=""):
    return old_primer.sequence.sequence == \
        tail + \
        getdna(old_primer.chromosome,old_primer.start_pos,old_primer.end_pos)

def check_right_primer(old_primer, tail=""):
    return old_primer.sequence.sequence == \
        tail + \
        rc(getdna(old_primer.chromosome,old_primer.start_pos,
            old_primer.end_pos))

def get_left_primer_company_tag(old_primer, tail=""):
    m = re.match("^{}([ACTG]){}$".format(
            tail,
            getdna(old_primer.chromosome,old_primer.start_pos,
                old_primer.end_pos),
        ),
        old_primer.sequence.sequence
    )
    if m:
        return m.groups()[0]

def get_right_primer_company_tag(old_primer, tail=""):
    m = re.match("^{}([ACTG]){}$".format(
            tail,
            rc(getdna(old_primer.chromosome,old_primer.start_pos,
                old_primer.end_pos)),
        ),
        old_primer.sequence.sequence
    )
    if m:
        return m.groups()[0]

def prepare_primer_dict(old_primer, ugs_model, apps, schema_editor):
    """
    Given an old_primer, and the name of the respective UGS model, prepares
    the common fields for the new primer.
    """
    db_alias = schema_editor.connection.alias
    UGS = apps.get_model("planning", ugs_model)
    slice = get_or_create_slice(old_primer, apps, schema_editor)
    ugs, c = UGS.objects.using(db_alias).get_or_create(
        slice=slice,
    )
    # NOTE: we rely on the user of this method to save, after updating
    # old_primer further.
    old_primer.new_ugs_id = ugs.id
    return {
        "name": old_primer.name,
        "ugs": ugs,
    }

def prepare_left_primer_dict(old_primer, apps, schema_editor):
    return prepare_primer_dict(old_primer, "UGSPlus", apps, schema_editor)

def prepare_right_primer_dict(old_primer, apps, schema_editor):
    return prepare_primer_dict(old_primer, "UGSMinus", apps, schema_editor)

def get_or_create_left_irac(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    IlluminaReadingAdaptor1 = apps.get_model("parts","IlluminaReadingAdaptor1")
    IlluminaReadingAdaptor1Cuts = apps.get_model("parts","IlluminaReadingAdaptor1Cuts")
    ira, c = IlluminaReadingAdaptor1.objects.using(db_alias).get_or_create(
        name="Illumina Standard Reading Adaptor1",
        _sequence=IRA1,
    )
    irac, c = IlluminaReadingAdaptor1Cuts.objects.using(db_alias).get_or_create(
        ira=ira,
        overlap_start=IRAC1_OVERLAP_START,
        overlap_end=IRAC1_OVERLAP_END,
    )
    return irac

def get_or_create_right_irac(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    IlluminaReadingAdaptor2 = apps.get_model("parts","IlluminaReadingAdaptor2")
    IlluminaReadingAdaptor2Cuts = apps.get_model("parts","IlluminaReadingAdaptor2Cuts")
    ira, c = IlluminaReadingAdaptor2.objects.using(db_alias).get_or_create(
        name="Illumina Standard Reading Adaptor2",
        _sequence=IRA2,
    )
    irac, c = IlluminaReadingAdaptor2Cuts.objects.using(db_alias).get_or_create(
        ira=ira,
        overlap_start=IRAC2_OVERLAP_START,
        overlap_end=IRAC2_OVERLAP_END,
    )
    return irac

@transaction.atomic
def convert_left_primer_tail(qs, apps, schema_editor):
    print
    print "Converting left primers with tails:"
    db_alias = schema_editor.connection.alias
    PCR1PlusPrimer = apps.get_model("synthesis","PCR1PlusPrimer")
    PCR1WithCompanyTagPlusPrimer = apps.get_model("synthesis","PCR1WithCompanyTagPlusPrimer")
    irac = get_or_create_left_irac(apps, schema_editor)
    for old_target in bar(qs):
        old_primer = old_target.primer
        d = prepare_left_primer_dict(old_primer, apps, schema_editor)
        d["irac"] = irac
        if check_left_primer(old_primer,tail=old_primer.tail.tail):
            primer = PCR1PlusPrimer.objects.using(db_alias).create(**d)
        else:
            tag = get_left_primer_company_tag(
                old_primer,tail=old_primer.tail.tail)
            if tag:
                d["tag"] = tag
                primer = PCR1WithCompanyTagPlusPrimer.objects.using(db_alias) \
                    .create(**d)
            else:
                raise IntegrityError("Bad left primer with tail: sequence and "
                    "descriptive fields don't match.")
        old_primer.new_primer_id = primer.id
        old_primer.new_primer_model = primer._meta.model_name
        old_primer.save()
        transfer_physical_locations(old_primer, primer, apps, schema_editor)

@transaction.atomic
def convert_right_primer_tail(qs, apps, schema_editor):
    print
    print "Converting right primers with tails:"
    db_alias = schema_editor.connection.alias
    PCR1MinusPrimer = apps.get_model("synthesis","PCR1MinusPrimer")
    PCR1WithCompanyTagMinusPrimer = apps.get_model("synthesis","PCR1WithCompanyTagMinusPrimer")
    irac = get_or_create_right_irac(apps, schema_editor)
    for old_target in bar(qs):
        old_primer = old_target.primer
        d = prepare_right_primer_dict(old_primer, apps, schema_editor)
        d["irac"] = irac
        if check_right_primer(old_primer,tail=old_primer.tail.tail):
            primer = PCR1MinusPrimer.objects.using(db_alias).create(**d)
        else:
            tag = get_right_primer_company_tag(
                old_primer,tail=old_primer.tail.tail)
            if tag:
                d["tag"] = tag
                primer = PCR1WithCompanyTagMinusPrimer.objects.using(db_alias) \
                    .create(**d)
            else:
                raise IntegrityError("Bad right primer with tail: sequence and "
                    "descriptive fields don't match.")
        old_primer.new_primer_id = primer.id
        old_primer.new_primer_model = primer._meta.model_name
        old_primer.save()
        transfer_physical_locations(old_primer, primer, apps, schema_editor)

@transaction.atomic
def convert_left_primer_notail(qs, apps, schema_editor):
    print
    print "Converting left primers without tails:"
    db_alias = schema_editor.connection.alias
    TargetedNoTailPlusPrimer = apps.get_model("synthesis","TargetedNoTailPlusPrimer")
    for old_target in bar(qs):
        old_primer = old_target.primer
        d = prepare_left_primer_dict(old_primer, apps, schema_editor)
        if check_left_primer(old_primer):
            primer = TargetedNoTailPlusPrimer.objects.using(db_alias).create(**d)
        else:
            raise IntegrityError("Bad left primer without tail: sequence and "
                "descriptive fields don't match.")
        old_primer.new_primer_id = primer.id
        old_primer.new_primer_model = primer._meta.model_name
        old_primer.save()
        transfer_physical_locations(old_primer, primer, apps, schema_editor)

@transaction.atomic
def convert_right_primer_notail(qs, apps, schema_editor):
    print
    print "Converting right primers without tails:"
    db_alias = schema_editor.connection.alias
    TargetedNoTailMinusPrimer = apps.get_model("synthesis","TargetedNoTailMinusPrimer")
    for old_target in bar(qs):
        old_primer = old_target.primer
        d = prepare_right_primer_dict(old_primer, apps, schema_editor)
        if check_right_primer(old_primer):
            primer = TargetedNoTailMinusPrimer.objects.using(db_alias).create(**d)
        else:
            raise IntegrityError("Bad right primer without tail: sequence and "
                "descriptive fields don't match.")
        old_primer.new_primer_id = primer.id
        old_primer.new_primer_model = primer._meta.model_name
        old_primer.save()
        transfer_physical_locations(old_primer, primer, apps, schema_editor)

@transaction.atomic
def convert_unknown_primer_notail(qs, apps, schema_editor):
    print
    print "Converting unknown primers with tails:"
    db_alias = schema_editor.connection.alias
    TargetedNoTailPlusPrimer = apps.get_model("synthesis","TargetedNoTailPlusPrimer")
    TargetedNoTailMinusPrimer = apps.get_model("synthesis","TargetedNoTailMinusPrimer")
    for old_target in bar(qs):
        old_primer = old_target.primer
        if check_left_primer(old_primer):
            d = prepare_left_primer_dict(old_primer, apps, schema_editor)
            primer = TargetedNoTailPlusPrimer.objects.using(db_alias).create(**d)
        elif check_right_primer(old_primer):
            d = prepare_right_primer_dict(old_primer, apps, schema_editor)
            primer = TargetedNoTailMinusPrimer.objects.using(db_alias).create(**d)
        else:
            raise IntegrityError("Bad primer without tail: sequence and "
                "descriptive fields don't match, neither as left nor as right.")
        old_primer.new_primer_id = primer.id
        old_primer.new_primer_model = primer._meta.model_name
        old_primer.save()
        transfer_physical_locations(old_primer, primer, apps, schema_editor)

def split_primers(qs, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    OldPrimerTail = apps.get_model("linapp", "PrimerTail")
    left_tail = OldPrimerTail.objects.using(db_alias).get(tail=LEFT_TAIL)
    right_tail = OldPrimerTail.objects.using(db_alias).get(tail=RIGHT_TAIL)
    queries = [
        (
            {
                "right_c": 0,
                "primer__tail": left_tail,
            },
            convert_left_primer_tail,
        ),
        (
            {
                "left_c": 0,
                "primer__tail": right_tail,
            },
            convert_right_primer_tail,
        ),
        (
            {
                "left_c__gt": 0,
                "primer__tail__isnull": True,
            },
            convert_left_primer_notail,
        ),
        (
            {
                "right_c__gt": 0,
                "primer__tail__isnull": True,
            },
            convert_right_primer_notail,
        ),
        (
            {
                "left_c": 0,
                "right_c": 0,
                "primer__tail__isnull": True,
            },
            convert_unknown_primer_notail,
        ),
    ]
    old_primer_targets = qs.select_related(
        "primer",
        "primer__sequence",
        "primer__tail",
        "chromosome",
    ).annotate(
            left_c=Count("primer__left_primer"),
            right_c=Count("primer__right_primer"),
    )
    try:
        check_partition(old_primer_targets,[q for q,f in queries])
    except NotCovering as e:
        raise IntegrityError("We have a wild Primer, example: {}".format(
            e.args[0][0].id))
    except NotMutuallyExclusive as e:
        raise IntegrityError("We have an ambiguous Primer - queries #{} and " \
            "#{} both match. Example: {}".format(e.args[0],e.args[1],
                e.args[2][0].id))
    for q,f in queries:
        f(old_primer_targets.filter(**q), apps, schema_editor)

