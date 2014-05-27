__author__ = 'ofirr'
from linapp.models import Sequence, Microsatellite
from math import floor
import hashlib
def trim_ms(ms):
    ms.repeat_number = floor(ms.repeat_number)
    seq = ms.referencevalue.sequence[:int(ms.repeat_unit_len*floor(ms.repeat_number))]
    try:
        sequence = Sequence.objects.get(hash=hashlib.md5(seq).hexdigest())
    except Sequence.DoesNotExist:
        sequence = Sequence(length=len(seq), sequence=seq, hash=hashlib.md5(seq).hexdigest())
        sequence.save()
    ms.end_pos = ms.start_pos + len(sequence.sequence) - 1
    assert ms.chromosome.getdna(ms.start_pos, ms.end_pos) == sequence.sequence
    ms.referencevalue = sequence
    ms.name = '{}_{}_{}'.format(ms.chromosome.name, ms.start_pos, ms.end_pos)
    ms.save()


def cut_ms(ms, k=1):
    ms.repeat_number = floor(ms.repeat_number)
    seq = ms.referencevalue.sequence[:int(ms.repeat_unit_len*(floor(ms.repeat_number)-k))]
    try:
        sequence = Sequence.objects.get(hash=hashlib.md5(seq).hexdigest())
    except Sequence.DoesNotExist:
        sequence = Sequence(length=len(seq), sequence=seq, hash=hashlib.md5(seq).hexdigest())
        sequence.save()
    ms.end_pos = ms.start_pos + len(sequence.sequence) - 1
    assert ms.chromosome.getdna(ms.start_pos, ms.end_pos) == sequence.sequence
    ms.referencevalue = sequence
    ms.name = '{}_{}_{}'.format(ms.chromosome.name, ms.start_pos, ms.end_pos)
    ms.save()



def resolve_overlapping_mss(te_queryset, max_trim=3):
    c = 0
    for te in te_queryset:
        for t in te.targets.all():
            c += 1
            ms = Microsatellite.objects.get(pk=t.pk)
            overlapping_mss = Microsatellite.objects.filter(chromosome=ms.chromosome).filter(start_pos__lte=ms.end_pos).filter(end_pos__gte=ms.end_pos).exclude(pk=ms.pk)
            if overlapping_mss:
                for ms2 in overlapping_mss:
                    if ms.start_pos == ms2.start_pos and ms.end_pos == ms2.end_pos:
                        ms.delete()
                        continue

                trim_ms(ms)
            overlapping_mss = Microsatellite.objects.filter(chromosome=ms.chromosome).filter(start_pos__lte=ms.end_pos).filter(end_pos__gte=ms.end_pos).exclude(pk=ms.pk)
            for i in range(max_trim):
                if overlapping_mss:
                    cut_ms(ms)
            overlapping_mss = Microsatellite.objects.filter(chromosome=ms.chromosome).filter(start_pos__lte=ms.end_pos).filter(end_pos__gte=ms.end_pos).exclude(pk=ms.pk)
            if overlapping_mss:
                print ms.repeat_unit_len, [oms.repeat_unit_len for oms in overlapping_mss]
                print ms, overlapping_mss
                print ms.repeat_number, [oms.repeat_number for oms in overlapping_mss]
                print ms.pk, [oms.pk for oms in overlapping_mss]
                print
            if c % 100 == 0:
                print '############'
                print c
                print '############'