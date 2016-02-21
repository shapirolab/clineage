import pytest
from models import Organ, Tissue, Individual, ExtractionEvent, Extraction, SamplingEvent, Cell, SampleComposition
from misc.test_models import human_taxa
from accounts.test_user import user
import datetime

@pytest.fixture()
def organ(db):
    o = Organ.objects.create(
        name="Blood",
    )
    return o


@pytest.fixture()
def tissue(db):
    t = Tissue.objects.create(
        name="Blood",
    )
    return t


@pytest.fixture()
def composition(db):
    sc = SampleComposition.objects.create(
        name="Single Cell"
    )
    return sc


@pytest.fixture()
def human_individual(human_taxa):
    i = Individual.objects.create(
        taxa=human_taxa,
        sex="M",
        name="Yossi",
    )
    return i


@pytest.fixture()
def human_extractionevent(human_individual, user):
    ee = ExtractionEvent.objects.create(
        individual=human_individual,
        name='human_extractionevent',
        date=datetime.datetime.now(),
        user_performed=user,
        user_documented=user,
    )
    return ee


@pytest.fixture()
def human_extraction(human_extractionevent, organ, tissue):
    e = Extraction.objects.create(
        extraction_event=human_extractionevent,
        name='human_extraction',
        organ=organ,
        tissue=tissue,
    )
    return e


@pytest.fixture()
def human_samplingevent(human_extraction):
    e = SamplingEvent.objects.create(
        extraction=human_extraction,
        name='human_samplingevent',
    )
    return e


@pytest.fixture()
def human_cell(human_individual, composition, human_samplingevent=None):
    c = Cell.objects.create(
        individual=human_individual,
        sampling=human_samplingevent,
        name='human_cell',
        composition=composition,
    )
    return c


@pytest.mark.django_db
def test_individual(human_individual):
    assert human_individual.name == "Yossi"


@pytest.mark.django_db
def test_cell_without_sampling(human_cell):
    assert human_cell.name == "human_cell"
    assert human_cell.sampling is None


@pytest.mark.django_db
def test_cell_with_sampling(human_samplingevent, composition):
    i = human_samplingevent.extraction.extraction_event.individual
    c = human_cell(i, composition, human_samplingevent)
    assert c.name == "human_cell"
    assert c.sampling.name == "human_samplingevent"
