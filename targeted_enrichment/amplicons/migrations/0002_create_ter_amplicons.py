# -*- coding: utf-8 -*-
# Generated manually, in form of Django 1.9.4 on 2016-05-02 19:56
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion
from misc.dna import DNA

TER_SELECT_RELATED = [
    "left_primer",
    "left_primer__ugs",
    "left_primer__ugs__slice",
    "right_primer",
    "right_primer__ugs",
    "right_primer__ugs__slice",
]


def get_target_amplicon_dict(ter, apps, schema_editor):
    db_alias = schema_editor.connection.alias
    DNASlice = apps.get_model('genomes', 'DNASlice')
    if ter.left_primer.ugs.slice.chromosome_id != \
        ter.right_primer.ugs.slice.chromosome_id:
        raise IntegrityError(
            "TER {} has ugss with slices from different " \
            "chromosomes.".format(ter))
    slice, c = DNASlice.objects.using(db_alias).get_or_create(
        chromosome_id=ter.left_primer.ugs.slice.chromosome_id,
        start_pos=ter.left_primer.ugs.slice.end_pos + 1,  # 1 based!
        end_pos=ter.right_primer.ugs.slice.start_pos - 1,  # 1 based!
    )
    return {
        "left_ugs": ter.left_primer.ugs,
        "right_ugs": ter.right_primer.ugs,
        "slice": slice,
    }


def populate_plain_targeted_ters(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    PCR1PrimerPairTER = apps.get_model('reagents', 'PCR1PrimerPairTER')
    PCR1PrimerPairTERDeprecated = apps.get_model('reagents', 'PCR1PrimerPairTERDeprecated')
    TargetedNoTailPrimerPairTER = apps.get_model('reagents', 'TargetedNoTailPrimerPairTER')
    PlainTargetedAmplicon = apps.get_model('amplicons', 'PlainTargetedAmplicon')
    for TER in [
        PCR1PrimerPairTER,
        PCR1PrimerPairTERDeprecated,
        TargetedNoTailPrimerPairTER,
    ]:
        for ter in TER.objects.using(db_alias).select_related(
                *TER_SELECT_RELATED):
            amplicon = PlainTargetedAmplicon.objects.using(db_alias).create(
                **get_target_amplicon_dict(ter, apps, schema_editor)
            )
            ter.amplicon = amplicon
            ter.save()


def populate_pcr1withcompanytagprimerpairter(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    PCR1WithCompanyTagPrimerPairTER = apps.get_model('reagents', 'PCR1WithCompanyTagPrimerPairTER')
    TargetedAmpliconWithCompanyTag = apps.get_model('amplicons', 'TargetedAmpliconWithCompanyTag')
    for ter in PCR1WithCompanyTagPrimerPairTER.objects.using(db_alias). \
        select_related(*TER_SELECT_RELATED):
        amplicon = TargetedAmpliconWithCompanyTag.objects.using(db_alias). \
            create(
            left_tag=ter.left_primer.tag,
            right_tag=DNA(ter.right_primer.tag).rev_comp().seq.decode("ascii"),
            **get_target_amplicon_dict(ter, apps, schema_editor)
        )
        ter.amplicon = amplicon
        ter.save()


def populate_shortpadlockfirstter(apps, schema_editor):
    db_alias = schema_editor.connection.alias
    OM6PadlockTER = apps.get_model('reagents', 'OM6PadlockTER')
    UMITargetedAmplicon = apps.get_model('amplicons', 'UMITargetedAmplicon')
    for ter in OM6PadlockTER.objects.using(db_alias). \
        select_related(
            "padlock",
            "padlock__left_ugs",
            "padlock__left_ugs__slice",
            "padlock__right_ugs",
            "padlock__right_ugs__slice",
        ):
        if ter.padlock.left_ugs.slice.chromosome_id != \
            ter.padlock.right_ugs.slice.chromosome_id:
            raise IntegrityError(
                "TER {} has ugss with slices from different " \
                "chromosomes.".format(ter))
        slice = DNASlice.objects.using(db_alias).create(
            chromosome_id=ter.padlock.left_ugs.slice.chromosome_id,
            start_pos=ter.padlock.left_ugs.slice.end_pos + 1,  # 1 based!
            end_pos=ter.padlock.right_ugs.slice.start_pos - 1,  # 1 based!
        )
        amplicon = UMITargetedAmplicon.objects.using(db_alias).create(
            umi_length=ter.padlock.umi_length,
            left_ugs=ter.padlock.left_ugs,
            right_ugs=ter.padlock.left_ugs,
            slice=slice,
        )
        ter.amplicon = amplicon
        ter.save()


def populate_ter_amplicons(apps, schema_editor):
    populate_plain_targeted_ters(apps, schema_editor)
    populate_pcr1withcompanytagprimerpairter(apps, schema_editor)
    populate_shortpadlockfirstter(apps, schema_editor)


class Migration(migrations.Migration):

    dependencies = [
        ('amplicons', '0001_initial'),
        ('reagents', '0004_auto_20160420_1723'),
    ]

    operations = [
        migrations.RunPython(
            code=populate_ter_amplicons,
            reverse_code=migrations.RunPython.noop,
        ),
    ]
