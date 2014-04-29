import re
from LinApp.models import *
from Utils.wells import *
from collections import defaultdict
from decimal import Decimal

def getTabDelimitedValues(f):
    content = f.read()
    f.close()
    lines = []
    for line in content.split('\n'):
        lines.append(line.split('\t'))
    return lines


class RelocRow(object):
    """
    Define table columns
    """
    name = 0
    reagent_type = 1
    volume = 2
    conc = 3
    plate_name = 4
    plate_type = 5
    well = 6
    dst_plate_name = 7
    dst_plate_type = 8
    dst_well = 9
    dst_volume = 10
    dest_conc = 11
    protocol = 12


MPX_LINE_LENGTH = 13


def get_rgts_map(tabDelimitedValues):
    rgts_map = defaultdict(list)
    for row in tabDelimitedValues[1:]:
        print 'proccessing row: %s' % row[RelocRow.name]
        if len(row) == MPX_LINE_LENGTH:
            try:
                plate = Plate.object.get(name=row[RelocRow.plate_name])
                location = SampleLocation.objects.get(plate=plate, well=row[RelocRow.well])
                reagent = location.reagent
                protocol = Protocol.objects.get(name=row[RelocRow.protocol]) if row[RelocRow.protocol] else None
                dst_plate_type = PlateType.objects.get(name=row[RelocRow.dst_plate_type])
            except Plate.DoesNotExist:
                print 'reference to non existing plate: %s' % row[RelocRow.plate_name]
                raise
            except SampleLocation.DoesNotExist:
                print 'reference to non existing location at plate: %s, location: %s' \
                      % (row[RelocRow.plate_name], row[RelocRow.well])
                raise
            except Protocol.DoesNotExist:
                print 'reference to non existing protocol: %s' % row[RelocRow.protocol]
                raise
            dst = {'dst_plate_name': row[RelocRow.dst_plate_name],
                   'dst_plate_type': dst_plate_type,
                   'dst_well': row[RelocRow.dst_well],
                   'dst_volume': row[RelocRow.dst_volume],
                   'dst_conc': row[RelocRow.dest_conc],
                   'protocol': protocol,}
            rgts_map[reagent] = dst
        else:
            print 'line does not match required length %d' % MPX_LINE_LENGTH
            raise
    return rgts_map


def relocate_rgts_map(rgts_map):
    for key in rgts_map.keys():
        plate, plate_created = Plate.objects.get_or_create(name=rgts_map[key]['dst_plate_name'],
                                                           type=rgts_map[key]['dst_plate_type'])
        location, location_created = SampleLocation.objects.get_or_create(plate=plate,
                                                                          well=rgts_map[key]['dst_well'])
        if not location_created:
            print 'Error: location %s %s already in use' % (plate, rgts_map[key]['dst_well'])
            raise
        location.reagent = key
        location.volume = Decimal(rgts_map[key]['dst_volume'])
        location.conc = Decimal(rgts_map[key]['dst_conc'])
        location.save()
