__author__ = 'ofirr'
import xlwt
import xlrd
from utils.wells import index2str
from linapp.models import SampleLocation, Assembly, Plate, PlateType


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def position_in_plates(targetenrichments_group,
                       stk_fw_plate,
                       stk_rv_plate,
                       pairs_plate,
                       plate_size=96,
                       pair_vol=50,
                       pair_conc=60,
                       stk_vol=100,
                       stk_conc=100):
    """
    Positions up to 96 primer pairs in stk plates and mixed pairs.
    returns three lists of (SampleLocation, created) tuples for pair, fw, rv
    """
    assert len(targetenrichments_group) <= plate_size
    fw_positions, rv_positions, pairs_positions = [], [], []
    for i, te in enumerate(targetenrichments_group):
        well = index2str(i+1)
        pair_location = SampleLocation.objects.get_or_create(
            plate=pairs_plate,
            well=well,
            defaults={'reagent': te,
                      'volume': pair_vol,
                      'concentration': pair_conc})
        pairs_positions.append(pair_location)
        stk_fw_location = SampleLocation.objects.get_or_create(
            plate=stk_fw_plate,
            well=well,
            defaults={'reagent': te.left,
                      'volume': stk_vol,
                      'concentration': stk_conc})
        fw_positions.append(stk_fw_location)
        stk_rv_location = SampleLocation.objects.get_or_create(
            plate=stk_rv_plate,
            well=well,
            defaults={'reagent': te.right,
                      'volume': stk_vol,
                      'concentration': stk_conc})
        rv_positions.append(stk_rv_location)
    return fw_positions, rv_positions, pairs_positions


def create_next_primers_plates(assembly):
    """
    This is a temporary function for generating the next primers plate
    """
    stk_primers_type = PlateType.objects.get(friendly='Primers STK')
    paired_primers_type = PlateType.objects.get(friendly='Primer pairs')
    mm9 = Assembly.objects.get(friendly_name='mm9')
    hg19 = Assembly.objects.get(friendly_name='hg19')
    if assembly == mm9:
        next_primers_plate_num = Plate.objects.filter(name__contains='United_MM_plate').count()+2
        name_united = 'United_MM_plate{}'.format(next_primers_plate_num)
        plate_united = Plate.objects.create(name=name_united, type=paired_primers_type)
        name_fw = 'MM_v2_Pl{}_Fw'.format(next_primers_plate_num)
        plate_fw = Plate.objects.create(name=name_fw, type=stk_primers_type)
        name_rev = 'MM_v2_Pl{}_Rev'.format(next_primers_plate_num)
        plate_rev = Plate.objects.create(name=name_rev, type=stk_primers_type)
    elif assembly == hg19:
        next_primers_plate_num = Plate.objects.filter(name__contains='United_hg19_Tails').count()+1
        name_united = 'United_hg19_Tails_plt{}'.format(next_primers_plate_num)
        plate_united = Plate.objects.create(name=name_united, type=paired_primers_type)
        name_fw = 'hg19_Tails_plt{}_Fw'.format(next_primers_plate_num)
        plate_fw = Plate.objects.create(name=name_fw, type=stk_primers_type)
        name_rev = 'hg19_Tails_plt{}_Rev'.format(next_primers_plate_num)
        plate_rev = Plate.objects.create(name=name_rev, type=stk_primers_type)
    else:
        print 'ERROR: unsupported assembly'
        raise
    return plate_united, plate_fw, plate_rev


def insertion_plates_to_db(create_primer_pairs, xls_name, plate_united, plate_fw, plate_rev, assembly='hg19', plate_size=96):
    primers_partition_list = list(chunks(range(1, len(create_primer_pairs)), plate_size))
    for plate in primers_partition_list:
        targetenrichments_group = create_primer_pairs.keys(plate)
        pairs_plate, stk_fw_plate, stk_rv_plate = create_next_primers_plates(assembly)
        fw_positions, rv_positions, pairs_positions = position_in_plates(targetenrichments_group,
                                                                       stk_fw_plate,
                                                                       stk_rv_plate,
                                                                       pairs_plate)
        create_primer_order_file_xls(targetenrichments_group, fw_positions, xls_name)
        create_primer_order_file_xls(targetenrichments_group, rv_positions, xls_name)
    return fw_positions, rv_positions, pairs_positions


def create_primer_order_file_xls(plate, positions, xls_name):

    workbook = xlwt.Workbook(xls_name)
    sheet = workbook.add_sheet(plate.name)
    for index, value in enumerate(positions):
        sheet.write(0, index, value)
    workbook.save(xls_name)
