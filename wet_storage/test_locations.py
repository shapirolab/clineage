import pytest
import datetime
from models import SampleLocation

from misc.test_models import human_taxa
from genomes.test_models import hg19_assembly, hg19_chromosome, \
    slice_28727_left, slice_28727_right, slice_28727_target_a, slice_28727_target_b,\
    slice_28734_left, slice_28734_right, slice_28734_target_a
from targeted_enrichment.planning.test_models import te_28727, te_28734, \
    ugs_28727_left, ugs_28727_right, ugs_28734_left, ugs_28734_right, \
    ms_28727_a, ms_28727_b, ms_28734_a
from targeted_enrichment.reagents.test_models import ter_28727, ter_28734
from primers.synthesis.test_models import primer_28727_left, primer_28727_right, \
    primer_28734_left, primer_28734_right
from primers.parts.test_models import illuminareadingadaptor1cuts, illuminareadingadaptor2cuts, \
    illuminareadingadaptor1, illuminareadingadaptor2
from test_models import plate_united, platetype_pairs, platecontext, plateplastica_pairs


@pytest.fixture()
def samplelocation_28727_united(plate_united, ter_28727):
    sl = SampleLocation.objects.create(
        plate=plate_united,
        well='A01',
        reagent=ter_28727,
        volume=60.0,
        concentration=50.0,
        timestamp=datetime.datetime(2014, 1, 15, 13, 51, 35),
    )
    return sl


@pytest.fixture()
def samplelocation_28734_united(plate_united, ter_28727):
    sl = SampleLocation.objects.create(
        plate=plate_united,
        well='A08',
        reagent=ter_28727,
        volume=60.0,
        concentration=50.0,
        timestamp=datetime.datetime(2014, 1, 15, 13, 51, 37),
    )
    return sl


@pytest.mark.django_db
def test_plate_united(plate_united):
    assert plate_united.name == "United_hg19_Tails_plt1"
    assert plate_united.type.plastic.columns == 12