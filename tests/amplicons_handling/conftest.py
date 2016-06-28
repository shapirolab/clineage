import pytest
from .parser_primer_file import get_case_from_columns, process_row
from .primer_design import primer3_design, bowtie2_design, sort_best_primers
from .primers_insertion import create_primers_in_db
from .positioning import insertion_plates_to_db, create_primer_order_file_xls
from django.core.management import setup_environ


@pytest.fixture()
def primer3_28734(ugs_28734):
    input_name = 'test_ugs_input.txt'
    output_name = 'test_ugs_output.txt'

    p = primer3_design(ugs_28734,
                       input_name,
                       output_name)
    # So our objects don't have "special" objects in fields
    return p


@pytest.fixture()
def choose_primer_pair(sam_file, target_primers):
    input_name = 'test_ugs_input.txt'
    output_name = 'test_ugs_output.txt'

    p = sort_best_primers(sam_file, target_primers)
    # So our objects don't have "special" objects in fields
    return p


@pytest.fixture()
def choose_primer_pair(sam_file, target_primers):
    input_name = 'test_ugs_input.txt'
    output_name = 'test_ugs_output.txt'

    p = sort_best_primers(sam_file, target_primers)
    # So our objects don't have "special" objects in fields
    return p


@pytest.fixture()
def choose_primer_pair(sam_file, target_primers):
    input_name = 'test_ugs_input.txt'
    output_name = 'test_ugs_output.txt'

    p = sort_best_primers(sam_file, target_primers)
    # So our objects don't have "special" objects in fields
    return p
