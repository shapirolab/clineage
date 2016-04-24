from utils.SequenceManipulations import complement
from frogress import bar
from linapp.migrations.mig_0004.utils import getdna
from targeted_enrichment.planning.models import Microsatellite

def ms_type_permutations(ms_type_str):
    for i in range(len(ms_type_str)):
        perm = ms_type_str[i:] + ms_type_str[:i]
        yield perm
        yield complement(perm)


def get_sequence(s):
    if s._sequence is not None:
        return s._sequence
    else:
        return getdna(s.chromosome, s.start_pos, s.end_pos)


def fix_mss_repeat_unit_type():
    for rut_dict in Microsatellite.objects.values('repeat_unit_type').distinct():
        rut = rut_dict['repeat_unit_type']
        repeat_type_variants = set(ms_type_permutations(rut))
        print('Migrating MSs of type {}'.format(rut))
        for ms in bar(Microsatellite.objects.filter(repeat_unit_type=rut)):
            repeat_unit_ref_seq = get_sequence(ms.slice)[:ms.repeat_unit_len].decode("utf-8")
            if repeat_unit_ref_seq not in repeat_type_variants:
                trut = sorted(ms_type_permutations(repeat_unit_ref_seq))[0]
                ms.repeat_unit_type = trut
                ms.save()
                print('Fixed {}'.format(ms))