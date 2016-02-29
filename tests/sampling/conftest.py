import pytest
import datetime

from sampling.models import Organ, Tissue, Individual, ExtractionEvent, Extraction, SamplingEvent, Cell, SampleComposition

from tests.misc.conftest import *
from tests.accounts.conftest import *


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
def human_cell_no_se(human_individual, composition):
    c = Cell.objects.create(
        individual=human_individual,
        name='human_cell_no_se',
        composition=composition,
    )
    return c


@pytest.fixture()
def human_cell_with_se(human_individual, composition, human_samplingevent):
    c = Cell.objects.create(
        individual=human_individual,
        sampling=human_samplingevent,
        name='human_cell_with_se',
        composition=composition,
    )
    return c


