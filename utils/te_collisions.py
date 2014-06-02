__author__ = 'ofirr'
from linapp.models import Chromosome, TargetEnrichment

def delete_overlapping_tes():
    tes_to_delete = []
    tes_conflict = []

    c = 0
    t = TargetEnrichment.objects.all().count()
    for chr in Chromosome.objects.all():
        for te1 in TargetEnrichment.objects.filter(chromosome=chr):
            for te2 in TargetEnrichment.objects.filter(chromosome=chr)\
                                               .filter(left__start_pos__lte=te1.right.end_pos)\
                                               .filter(right__end_pos__gte=te1.left.start_pos)\
                                               .exclude(pk=te1.pk):
                if te1.physical_locations.all() and not te2.physical_locations.all():
                    te2.delete()
                    print 'deleted', te2
                    continue
                if te2.physical_locations.all() and not te1.physical_locations.all():
                    te1.delete()
                    print 'deleted', te1
                    continue
                if te2.physical_locations.all() and te1.physical_locations.all():
                    if te1.type != te2.type:
                        continue
                    tes_conflict.append((te1, te2))
                    continue
                if te1.left.start_pos < te2.left.start_pos and te1.right.end_pos > te2.right.end_pos:
                    te2.delete()
                    print 'deleted', te2
                    continue
                if te2.left.start_pos < te1.left.start_pos and te2.right.end_pos > te1.right.end_pos:
                    te1.delete()
                    print 'deleted', te1
                    continue
                tes_conflict.append((te1, te2))
            c += 1
            if c % 100 == 0:
                print len(tes_to_delete), len(tes_conflict), c, float(c)/t*100
    return tes_conflict