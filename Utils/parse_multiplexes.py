
import re
from LinApp.models import *
import hashlib
from collections import defaultdict
from django.db.models import Max, Count
from Utils.wells import *
from Utils.distribute_multiplex import *


def test_multiplex():
    f = open(r'C:\Users\ofirr\Desktop\human_mpx_23.10.2013.txt', 'rb')
    tabDeleimitedValues = getTabDelimitedValues(f)
    map = get_mpxs_map(tabDeleimitedValues)
    objs = create_mpx_objects_from_mpxs_map(map)
    locations = position_multiplexes(objs)
    res = []
    for obj in objs:
        target = obj.physical_locations.all()[0]
        targetEnrichments = obj.primers.all()
        sources = []
        for t in targetEnrichments:
            try:
                sources.append(t.physical_locations.all()[0])
            except Exception as e:
                print 'Error: exception in target %s:\r\n\t%s' % (str(t), str(e))
        for source in sources:
            tup = ((source.plate.name, source.well), (target.plate.name, target.well))
            res.append(tup)
    print res
    distribute_multiplex(res)


def getTabDelimitedValues(f):
    content = f.read()
    f.close()
    lines = []
    for line in content.split('\n'):
        lines.append(line.split('\t'))
    return lines


class MpxRow(object):
    """
    Define table columns
    """
    name = 0
    ignor1 = 1
    ignor2 = 2
    ignor3 = 3
    start_index = 4
    end_index = 5
    ignor6 = 6
    ignor7 = 7
    ignor8 = 8
    ignor9 = 9
    ignor10 = 10
    ignor11 = 11
    ignor12 = 12
    ignor13 = 13
    mpx = 14


MPX_LINE_LENGTH = 15


def get_mpxs_map(tabDelimitedValues):
    mpxs_map = defaultdict(list)
    for row in tabDelimitedValues[1:]:
        print 'proccessing row: %s' % row[MpxRow.name]
        if len(row) == MPX_LINE_LENGTH:
            assert re.match('^[0-9]+$', row[MpxRow.start_index])
            assert re.match('^[0-9]+$', row[MpxRow.end_index])
            try:
                target = Target.objects.get(name=row[MpxRow.name],
                                            start_pos=row[MpxRow.start_index],
                                            end_pos=row[MpxRow.end_index])
            except Target.DoesNotExist:
                print 'reference to non existing MS: %s' % row[MpxRow.name]
                raise
            if row[MpxRow.mpx]:
                mpxs_map[row[MpxRow.mpx]].append(target)
    return mpxs_map


def create_mpx_objects_from_mpxs_map(mpxs_map):
    mpx_objects = []
    keys = list(mpxs_map.keys())
    keys.sort(key=lambda s: int(s))
    for mpx in keys:
        if PrimersMultiplex.objects.all():
            mpx_name = str(PrimersMultiplex.objects.all().aggregate(Max('pk')).get('pk__max')+1) + '_' + str(len(mpxs_map[mpx]))
        else:
            mpx_name = '1_' + str(len(mpxs_map[mpx]))
        mpx_object = PrimersMultiplex.objects.create(name=mpx_name)
        enrichments = []
        for target in mpxs_map[mpx]:
            dedicated_enrichments = target.targetenrichment_set.annotate(targets_count=Count('targets')).filter(targets_count=1)
            if len(dedicated_enrichments) == 1:
                enrichments.append(dedicated_enrichments[0])
            elif len(dedicated_enrichments) > 1:
                name_filtered = dedicated_enrichments\
                    .filter(left__name=target.name + '_fwd')\
                    .filter(right__name=target.name + '_rev')
                if len(name_filtered) == 1:
                    enrichments.append(name_filtered[0])
                else:
                    print 'Warning: multiple dedicated enrichments for target %s' % str(target)
                    enrichments.append(dedicated_enrichments[0])
            else:
                nonspecific_enrichments = target.targetenrichment_set.all()
                if len(nonspecific_enrichments) == 1:
                    print 'Warning: non specific enrichment for target %s' % str(target)
                    enrichments.append(nonspecific_enrichments[0])
                elif len(dedicated_enrichments) > 1:
                    print 'Warning: multiple non specific enrichments for target %s' % str(target)
                    enrichments.append(nonspecific_enrichments[0])
                else:
                    print 'Error: no enrichments found for target %s' % str(target)
                    raise
        mpx_object.primers = enrichments
        mpx_object.save()
        mpx_objects.append(mpx_object)
    return mpx_objects


def position_multiplexes(mpx_objects):
    mpx_plate_type = PlateType.objects.get(friendly='MPXs')
    next_id = Plate.objects.filter(type=mpx_plate_type).count() + 1
    mpx_plate = Plate.objects.create(type=mpx_plate_type,
                                 name='MM_MPX_w/o_LB_P'+str(next_id),
                                 lastusedwell='A01',
                                 barcode='',
                                 state='')
    locations = []
    for mpx in mpx_objects:
        location = SampleLocation.objects.create(plate=mpx_plate,
                                      well=index2str(str2index(mpx_plate.lastusedwell)),
                                      reagent=mpx,
                                      )
        mpx_plate.lastusedwell = index2str(str2index(mpx_plate.lastusedwell)+1)
        mpx_plate.save()
        locations.append(location)
    return locations

