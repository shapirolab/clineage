import pytest
from amplicons_handling.primer_design import primer3_design, sort_best_primers
from amplicons_handling.primers_insertion import create_target_enrichment_in_db
from amplicons_handling.positioning import insertion_OM_to_db
from tests.genomes.conftest import *
from tests.primers.parts.conftest import *
from tests.targeted_enrichment.planning.conftest import *
from targeted_enrichment.planning.models import Target
from primers.synthesis.models import PadlockPrepCommonPrimers


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
    chosen_target_primers, discarded_targets = sort_best_primers(primer3_output,
                                                                 size_filter=True, deltag_filter=False, delta_min=-85)

    return chosen_target_primers


@pytest.fixture()
def tate_tuple_sample(target_enrichment_sample):
    ta, te = create_target_enrichment_in_db(target_enrichment_sample)
    tate_tuple = (ta, te)

    return [tate_tuple]


@pytest.fixture()
def padlock_prep_primers(padlockamplificationplusprimerpart1, padlockamplificationplusprimerpart2,
                        padlockamplificationminusprimerpart1, padlockamplificationminusprimerpart2,
                         mly1_restriction_site_sample):
    padlock_prep_primers = PadlockPrepCommonPrimers.objects.create(
        name='sample_padlock',
        left_amp_primer_part1=padlockamplificationplusprimerpart1,
        left_amp_primer_part2=padlockamplificationplusprimerpart2,
        right_amp_primer_part1=padlockamplificationminusprimerpart1,
        right_amp_primer_part2=padlockamplificationminusprimerpart2,
        restriction_enzyme=mly1_restriction_site_sample,
    )
    # So our objects don't have "special" objects in fields
    padlock_prep_primers = PadlockPrepCommonPrimers.objects.get(pk=padlock_prep_primers.pk)
    return padlock_prep_primers


@pytest.fixture()
def OMmix_sample(tate_tuple_sample, illuminareadingadaptor1fortail, illuminareadingadaptor2fortail, padlock_prep_primers):
    OMmix = insertion_OM_to_db(tate_tuple_sample, 'sample_panel', illuminareadingadaptor1fortail, illuminareadingadaptor2fortail, 'sample_padlock')

    return OMmix
