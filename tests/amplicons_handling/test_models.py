import pytest
import xlrd
from misc.utils import get_unique_path, unlink
from amplicons_handling.primer_design import primer3_design, create_amplicon_for_primer3
from amplicons_handling.primers_insertion import create_target_enrichment_in_db
from amplicons_handling.positioning import insertion_OM_to_db, create_primer_order_file_xls


@pytest.mark.django_db
def test_create_amplicon_for_primer3(slice_28727_amplicon):
    slice = create_amplicon_for_primer3(slice_28727_amplicon, 1000)

    assert slice.start_pos == 81316117-1000
    assert slice.end_pos == 81316242+1000
    assert slice.chromosome == slice_28727_amplicon.chromosome


@pytest.mark.django_db
def test_get_primer3_output(sample_target):
    target_primers = primer3_design(sample_target)
    assert target_primers['SEQUENCE_ID'] == sample_target.id
    assert 130 <= target_primers['PRIMER_PAIR_1_PRODUCT_SIZE'] <= 300



@pytest.mark.django_db
def test_create_target_enrichment_in_db(target_enrichment_sample):
    ta_sample, te_sample = create_target_enrichment_in_db(target_enrichment_sample)

    assert target_enrichment_sample['RIGHT_START'][0] - \
                target_enrichment_sample['LEFT_START'][0] == \
           (te_sample.right.slice.end_pos - te_sample.left.slice.start_pos)
    assert te_sample.left.sequence.seq.decode('utf-8') == target_enrichment_sample['LEFT']
    assert te_sample.right.sequence.seq.decode('utf-8') == target_enrichment_sample['RIGHT']
    assert ta_sample.slice.start_pos == 81316078  # TODO: revise


@pytest.mark.django_db
def test_insertion_OM_to_db(tate_tuple_sample, illuminareadingadaptor1fortail, illuminareadingadaptor2fortail, padlock_prep_primers):
    sample_mix = insertion_OM_to_db(tate_tuple_sample, 'sample_panel', illuminareadingadaptor1fortail, illuminareadingadaptor2fortail, 'sample_padlock')
    om6_prep = sample_mix[0]
    ta, te = tate_tuple_sample[0]
    assert te.left == om6_prep.padlock.left_ugs
    assert om6_prep.primers.restriction_enzyme.name == 'Mly1'


@pytest.mark.django_db
def test_order_file_xls(OMmix_sample):

    with unlink(get_unique_path(ext='xls')) as sample_sheet:
        create_primer_order_file_xls(OMmix_sample, sample_sheet)

        workbook = xlrd.open_workbook(sample_sheet)
        worksheet = workbook.sheet_by_name('OM7')

        row_0 = worksheet.row(0)
        row_1 = worksheet.row(1)

    assert row_0[0].value == 'TEID'

