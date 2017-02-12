import pytest
import datetime

from sampling.models import Organ, Tissue, Individual, ExtractionEvent, Extraction, SamplingEvent, Cell, SampleComposition

from tests.misc.conftest import *
from tests.accounts.conftest import *


@pytest.fixture()
def organ(transactional_db):
    o = Organ.objects.create(
        name="Blood",
    )
    # So our objects don't have "special" objects in fields
    o = Organ.objects.get(pk=o.pk)
    return o


@pytest.fixture()
def tissue(transactional_db):
    t = Tissue.objects.create(
        name="Blood",
    )
    # So our objects don't have "special" objects in fields
    t = Tissue.objects.get(pk=t.pk)
    return t


@pytest.fixture()
def composition(transactional_db):
    sc = SampleComposition.objects.create(
        name="Single Cell"
    )
    # So our objects don't have "special" objects in fields
    sc = SampleComposition.objects.get(pk=sc.pk)
    return sc


@pytest.fixture()
def human_individual(human_taxa, user):
    i = Individual.objects.create(
        taxa=human_taxa,
        sex="M",
        name="Yossi",
        partner=user,
    )
    # So our objects don't have "special" objects in fields
    i = Individual.objects.get(pk=i.pk)
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
    # So our objects don't have "special" objects in fields
    ee = ExtractionEvent.objects.get(pk=ee.pk)
    return ee


@pytest.fixture()
def human_extraction(human_extractionevent, organ, tissue):
    e = Extraction.objects.create(
        extraction_event=human_extractionevent,
        name='human_extraction',
        organ=organ,
        tissue=tissue,
    )
    # So our objects don't have "special" objects in fields
    e = Extraction.objects.get(pk=e.pk)
    return e


@pytest.fixture()
def human_samplingevent(human_extraction):
    e = SamplingEvent.objects.create(
        extraction=human_extraction,
        name='human_samplingevent',
    )
    # So our objects don't have "special" objects in fields
    e = SamplingEvent.objects.get(pk=e.pk)
    return e


@pytest.fixture()
def human_cell_no_se(human_individual, composition):
    c = Cell.objects.create(
        individual=human_individual,
        name='human_cell_no_se',
        composition=composition,
    )
    # So our objects don't have "special" objects in fields
    c = Cell.objects.get(pk=c.pk)
    return c


@pytest.fixture()
def human_cell_with_se(human_individual, composition, human_samplingevent):
    c = Cell.objects.create(
        individual=human_individual,
        sampling=human_samplingevent,
        name='human_cell_with_se',
        composition=composition,
    )
    # So our objects don't have "special" objects in fields
    c = Cell.objects.get(pk=c.pk)
    return c
