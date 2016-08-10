__author__ = 'ofirr'
import xlwt
import xlrd
from wet_storage.models import SampleLocation
from lib_prep.multiplexes.models import OM6Oligomix
from targeted_enrichment.reagents.models import OM6PadlockTER
from primers.synthesis.models import OM6Padlock
from primers.parts.models import IlluminaReadingAdaptor1ForTail, IlluminaReadingAdaptor2ForTail


def insertion_OM_to_db(tate_tuple, panel_name):
    """
    creating list of ters for OM6 panel
    :param tate_tuple: UMITargetedAmplicon, TargetEnrichment list of tuples
    :panel_name: name for the panel
    :return:
    """

    ira1ft = IlluminaReadingAdaptor1ForTail.objects.filter(id=2)
    ira2ft = IlluminaReadingAdaptor2ForTail.objects.filter(id=2)
    OMmix = OM6Oligomix.objects.create(name=panel_name)
    ters = []
    for tate in tate_tuple:
        ta, te = tate
        om6_padlock = OM6Padlock.objects.get_or_create(
            left_ugs=te.left,
            right_ugs=te.right,
            ira1ft=ira1ft,
            ira2ft=ira2ft,
            umi_length=3,

        )
        ter = OM6PadlockTER.objects.get_or_create(
            te=te,
            amplicon=ta,
            padlock=om6_padlock,
        )
        ters.append(ter)
    OMmix.ters = ters

    return OMmix

def create_next_primers_plates(assembly):
    """
    This is a temporary function for generating the next primers plate
    """
    stk_primers_type = PlateType.objects.get(friendly='Primers STK')
    paired_primers_type = PlateType.objects.get(friendly='Primer pairs')
    #mm9 = Assembly.objects.get(friendly_name='mm9')
    #hg19 = Assembly.objects.get(friendly_name='hg19')
    if assembly == 'mm9':
        next_primers_plate_num = Plate.objects.filter(name__contains='United_MM_plate').count()+2
        name_united = 'United_MM_plate{}'.format(next_primers_plate_num)
        plate_united = Plate.objects.create(name=name_united, type=paired_primers_type)
        name_fw = 'MM_v2_Pl{}_Fw'.format(next_primers_plate_num)
        plate_fw = Plate.objects.create(name=name_fw, type=stk_primers_type)
        name_rev = 'MM_v2_Pl{}_Rev'.format(next_primers_plate_num)
        plate_rev = Plate.objects.create(name=name_rev, type=stk_primers_type)
    elif assembly == 'hg19':
        next_primers_plate_num = Plate.objects.filter(name__contains='United_hg19_Tails').count()+1
        name_united = 'United_hg19_Tails_plt{}'.format(next_primers_plate_num)
        plate_united = Plate.objects.create(name=name_united, type=paired_primers_type)
        name_fw = 'hg19_Tails_plt{}_Fw'.format(next_primers_plate_num)
        plate_fw = Plate.objects.create(name=name_fw, type=stk_primers_type)
        name_rev = 'hg19_Tails_plt{}_Rev'.format(next_primers_plate_num)
        plate_rev = Plate.objects.create(name=name_rev, type=stk_primers_type)
    else:
        print('ERROR: unsupported assembly')
        raise
    return plate_united, plate_fw, plate_rev


def insertion_plates_to_db(te_list, assembly='hg19', plate_size=96):
    pairs_plates, stk_fw_plates, stk_rv_plates = [], [], []
    for plate_te in chunks(te_list, plate_size):
        pairs_plate, stk_fw_plate, stk_rv_plate = create_next_primers_plates(assembly)
        print(pairs_plate)
        print(stk_fw_plate)
        print(stk_rv_plate)
        pairs_plates.append(pairs_plate)
        stk_fw_plates.append(stk_fw_plate)
        stk_rv_plates.append(stk_rv_plate)
        fw_positions, rv_positions, pairs_positions = position_in_plates(plate_te,
                                                                       stk_fw_plate,
                                                                       stk_rv_plate,
                                                                       pairs_plate)
    return pairs_plates, stk_fw_plates, stk_rv_plates


def create_primer_order_file_xls(stk_fw_plates, stk_rv_plates, xls_name):
    workbook = xlwt.Workbook()
    for fw_plate, rv_plate in zip(stk_fw_plates, stk_rv_plates):
        for plate in [fw_plate, rv_plate]:
            sheet = workbook.add_sheet(plate.name)
            sheet.write(0, 0, 'WellPosition')
            sheet.write(0, 1, 'Name')
            sheet.write(0, 2, 'Sequence')
            sheet.write(0, 3, 'Notes')
            for index, sl in enumerate(SampleLocation.objects.filter(plate=plate)):
                well = sl.well
                primer = sl.reagent
                name = primer.name
                primer_sequence = primer.sequence.sequence
                sheet.write(index+1, 0,  well)
                sheet.write(index+1, 1,  name)
                sheet.write(index+1, 2,  primer_sequence)
    workbook.save(xls_name)
