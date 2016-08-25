__author__ = 'ofirr'
import xlwt
import xlrd
from wet_storage.models import SampleLocation
from lib_prep.multiplexes.models import OM6Oligomix
from targeted_enrichment.reagents.models import OM6PadlockTER
from primers.synthesis.models import OM6Padlock, OM6Prep, PadlockPrepCommonPrimers


def insertion_OM_to_db(tate_tuple, panel_name, ira1ft, ira2ft):
    """
    creating list of ters for OM6 panel
    :param tate_tuple: UMITargetedAmplicon, TargetEnrichment list of tuples
    :panel_name: name for the panel
    :return:
    """

    OMmix = OM6Oligomix.objects.create(name=panel_name)
    ters = []
    prep_list = []
    for tate in tate_tuple:
        ta, te = tate
        om6_padlock, c = OM6Padlock.objects.get_or_create(
            left_ugs=te.left,
            right_ugs=te.right,
            ira1ft=ira1ft,
            ira2ft=ira2ft,
            umi_length=3,

        )
        om6_primers, c = PadlockPrepCommonPrimers.objects.get(name='OM6 Standard Adaptors')
        om6_prep, c = OM6Prep.objects.get_or_create(
            name=te.id,
            padlock=om6_padlock,
            primers=om6_primers,
        )
        prep_list.append(om6_prep)
        ter, c = OM6PadlockTER.objects.get_or_create(
            te=te,
            amplicon=ta,
            padlock=om6_padlock,
        )
        ters.append(ter)
    OMmix.ters = ters

    return prep_list


def create_primer_order_file_xls(prep_list, xls_name):
    workbook = xlwt.Workbook()

    sheet = workbook.add_sheet('OM7')
    sheet.write(0, 0, 'TEID')
    sheet.write(0, 1, 'Sequence')

    for index, prep in enumerate(prep_list):

        name = prep.name
        primer_sequence = prep.padlock.sequence.seq.decode('utf-8')
        sheet.write(index+1, 0,  name)
        sheet.write(index+1, 1,  primer_sequence)
    workbook.save(xls_name)


