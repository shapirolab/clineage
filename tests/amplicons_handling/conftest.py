import pytest
from amplicons_handling.primer_design import primer3_design, sort_best_primers
from amplicons_handling.primers_insertion import create_target_enrichment_in_db
from amplicons_handling.positioning import insertion_OM_to_db
from tests.genomes.conftest import *
from targeted_enrichment.planning.models import Target


@pytest.fixture()
def sample_target(slice_28727_amplicon):
    target_28727 = Target.objects.create(
        name='28727_test',
        slice=slice_28727_amplicon,
    )
    # So our objects don't have "special" objects in fields
    target_28727 = Target.objects.get(name=target_28727.name)
    return target_28727


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
