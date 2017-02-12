import pytest

from wet_storage.models import PlateContext, PlatePlastica, PlateType, Plate, SampleLocation
import datetime

@pytest.fixture()
def platecontext(transactional_db):
    pc = PlateContext.objects.create(
        description='Storage'
    )
    # So our objects don't have "special" objects in fields
    pc = PlateContext.objects.get(pk=pc.pk)
    return pc


@pytest.fixture()
def plateplastica_pairs(transactional_db):
    pp = PlatePlastica.objects.create(
        description='Tetrad',
        rows=8,
        columns=12,
    )
    # So our objects don't have "special" objects in fields
    pp = PlatePlastica.objects.get(pk=pp.pk)
    return pp


@pytest.fixture()
def plateplastica_stk(transactional_db):
    pp = PlatePlastica.objects.create(
        description='Deep-Well Square',
        rows=8,
        columns=12,
    )
    # So our objects don't have "special" objects in fields
    pp = PlatePlastica.objects.get(pk=pp.pk)
    return pp


@pytest.fixture()
def platetype_pairs(platecontext, plateplastica_pairs):
    pt = PlateType.objects.create(
        friendly='Primer pairs',
        context=platecontext,
        plastic=plateplastica_pairs,
    )
    # So our objects don't have "special" objects in fields
    pt = PlateType.objects.get(pk=pt.pk)
    return pt


@pytest.fixture()
def platetype_stk(platecontext, plateplastica_stk):
    pt = PlateType.objects.create(
        friendly='Primers STK',
        context=platecontext,
        plastic=plateplastica_stk,
    )
    # So our objects don't have "special" objects in fields
    pt = PlateType.objects.get(pk=pt.pk)
    return pt


@pytest.fixture()
def plate_united(platetype_pairs):
    p = Plate.objects.create(
        type=platetype_pairs,
        name='United_hg19_Tails_plt1',
    )
    # So our objects don't have "special" objects in fields
    p = Plate.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def plate_fwd(platetype_stk):
    p = Plate.objects.create(
        type=platetype_stk,
        name='hg19_Tails_plt1_Fw',
    )
    # So our objects don't have "special" objects in fields
    p = Plate.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def plate_rev(platetype_stk):
    p = Plate.objects.create(
        type=platetype_stk,
        name='hg19_Tails_plt1_Rev',
    )
    # So our objects don't have "special" objects in fields
    p = Plate.objects.get(pk=p.pk)
    return p


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
    # So our objects don't have "special" objects in fields
    sl = SampleLocation.objects.get(pk=sl.pk)
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
    # So our objects don't have "special" objects in fields
    sl = SampleLocation.objects.get(pk=sl.pk)
    return sl


@pytest.fixture()
def platetype_cells(platecontext, plateplastica_pairs):
    pt = PlateType.objects.create(
        friendly='Primer pairs',
        context=platecontext,
        plastic=plateplastica_pairs,
    )
    # So our objects don't have "special" objects in fields
    pt = PlateType.objects.get(pk=pt.pk)
    return pt


@pytest.fixture()
def cells_plate(platetype_cells):
    p = Plate.objects.create(
        type=platetype_cells,
        name='Human cells plate',
    )
    # So our objects don't have "special" objects in fields
    p = Plate.objects.get(pk=p.pk)
    return p


@pytest.fixture()
def samplelocation_human_cell_with_se(cells_plate, amplifiedcontent):
    sl = SampleLocation.objects.create(
        plate=cells_plate,
        well='A01',
        reagent=amplifiedcontent,
        timestamp=datetime.datetime(2015, 1, 1, 13, 51, 35),
    )
    # So our objects don't have "special" objects in fields
    sl = SampleLocation.objects.get(pk=sl.pk)
    return sl

@pytest.fixture()
def samplelocation_human_cell_no_se(cells_plate, amplifiedcontent_no_se):
    sl = SampleLocation.objects.create(
        plate=cells_plate,
        well='A01',
        reagent=amplifiedcontent_no_se,
        timestamp=datetime.datetime(2015, 1, 1, 13, 51, 35),
    )
    # So our objects don't have "special" objects in fields
    sl = SampleLocation.objects.get(pk=sl.pk)
    return sl
