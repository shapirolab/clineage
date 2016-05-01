from misc.dna import DNA
from frogress import bar
from targeted_enrichment.planning.models import Microsatellite

def ms_type_permutations(ms_type_str):
    for i in range(len(ms_type_str)):
        perm = ms_type_str[i:] + ms_type_str[:i]
        yield perm
        yield perm.rev_comp()


def fix_mss_repeat_unit_type():
    rut_canon = {}
    for ms in bar(Microsatellite.objects.all()):
        repeat_unit_ref_seq = ms.slice.sequence[:ms.repeat_unit_len]
        if repeat_unit_ref_seq not in rut_canon:
            repeat_type_variants = set(ms_type_permutations(repeat_unit_ref_seq))
            canon = min(repeat_type_variants, key=lambda x: x.seq)
            for in_rut in repeat_type_variants:
                rut_canon[in_rut] = canon
        canon = rut_canon[repeat_unit_ref_seq].seq.decode("ascii")
        if canon != ms.repeat_unit_type:
            old_ms = "{}".format(ms)
            ms.repeat_unit_type = canon
            ms.save()
            print('Fixed {} to {}'.format(old_ms, ms))
