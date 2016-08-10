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
def target_list_sample(sample_target):
    target_list = [sample_target]
    return target_list


@pytest.fixture()
def primer3_output(sample_target):
    target_primers = primer3_design(sample_target, primer_num_rerun=2)

    return target_primers


@pytest.fixture()
def target_enrichment_sample(primer3_output):
    chosen_target_primers, discarded_targets = sort_best_primers(primer3_output, delta_min=-85)

    return chosen_target_primers


@pytest.fixture()
def tate_tuple_sample(target_enrichment_sample):
    ta, te = create_target_enrichment_in_db(target_enrichment_sample)
    tate_tuple = (ta, te)

    return [tate_tuple]


@pytest.fixture()
def OMmix_sample(tate_tuple_sample):
    OMmix = insertion_OM_to_db(tate_tuple_sample, 'sample_panel')

    return OMmix
