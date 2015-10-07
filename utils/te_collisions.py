__author__ = 'ofirr'
from genomes.models import TargetEnrichment

def delete_overlapping_tes(te_queryset):
    tes_conflict = []
    tes_bug = []
    c = 0
    t = te_queryset.count()
    for te1 in te_queryset:
        for te2 in TargetEnrichment.objects.filter(chromosome=te1.chromosome)\
                                           .filter(left__start_pos__lte=te1.right.end_pos)\
                                           .filter(right__end_pos__gte=te1.left.start_pos)\
                                           .exclude(pk=te1.pk):
            try:
                if te1.physical_locations.all() and not te2.physical_locations.all():
                    print 'deleted', te2
                    te2.delete()
                    continue
                if te2.physical_locations.all() and not te1.physical_locations.all():
                    print 'deleted', te1
                    te1.delete()
                    continue
                if te2.physical_locations.all() and te1.physical_locations.all():
                    if te1.type != te2.type:
                        continue
                    tes_conflict.append((te1, te2))
                    continue
                if te1.left.start_pos < te2.left.start_pos and te1.right.end_pos > te2.right.end_pos:
                    print 'deleted', te2
                    te2.delete()
                    continue
                if te2.left.start_pos < te1.left.start_pos and te2.right.end_pos > te1.right.end_pos:
                    print 'deleted', te1
                    te1.delete()
                    continue
                tes_conflict.append((te1, te2))
            except AssertionError:
                tes_bug.append((te1,te2))
        c += 1
        if c % 100 == 0:
            print len(tes_conflict), len(tes_bug), c, float(c)/t*100
    return tes_conflict, tes_bug