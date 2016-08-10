import pytest
import xlrd
from amplicons_handling.primer_design import primer3_design, create_amplicon_for_primer3
from amplicons_handling.primers_insertion import create_target_enrichment_in_db
from amplicons_handling.positioning import insertion_OM_to_db, create_primer_order_file_xls


@pytest.mark.django_db
def test_create_amplicon_for_primer3(slice_28727_amplicon):
    slice = create_amplicon_for_primer3(slice_28727_amplicon, 1000)

    assert slice.start_pos == 81316117-1000
    assert slice.end_pos == 81316242+1000
    assert slice.chromosome == slice_28727_amplicon.chromosome
