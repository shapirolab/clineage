__author__ = 'ofirr'
from linapp.models import Chromosome, TargetEnrichment

def overlapping_tes():
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
                    tes_to_delete.append(te2)
                    #print 'deleted', te2
                    continue
                if te2.physical_locations.all() and not te1.physical_locations.all():
                    tes_to_delete.append(te1)
                    #print 'deleted', te
                    continue
                if te2.physical_locations.all() and te1.physical_locations.all():
                    print 'both exist', te1.physical_locations.all(), te2.physical_locations.all()
                    if te1.left.start_pos < te2.left.start_pos and te1.right.end_pos > te2.right.end_pos:
                        print 'te includes te2, te2 deleted', te1, te2
                        tes_conflict.append((te1, te2))
                        continue
                    if te2.left.start_pos < te1.left.start_pos and te2.right.end_pos > te1.right.end_pos:
                        print 'te2 includes te, te deleted', te1, te2
                        tes_conflict.append((te1, te2))
                        continue
                    tes_conflict.append((te1, te2))
                    continue
                if te1.left.start_pos < te2.left.start_pos and te1.right.end_pos > te2.right.end_pos:
                    tes_to_delete.append(te2)
                    continue
                if te2.left.start_pos < te1.left.start_pos and te2.right.end_pos > te1.right.end_pos:
                    tes_to_delete.append(te1)
                    continue
                print "#################shouldn't happen#################"
            c += 1
            if c % 100 == 0:
                print c, c/t*100
    return tes_to_delete, tes_conflict