import pytest
from models import PlateContext, PlatePlastica, PlateType, Plate


@pytest.fixture()
def platecontext(db):
    pc = PlateContext.objects.create(
        description='Storage'
    )
    return pc


@pytest.fixture()
def plateplastica_pairs(db):
    pp = PlatePlastica.objects.create(
        description='Tetrad',
        rows=8,
        columns=12,
    )
    return pp


@pytest.fixture()
def plateplastica_stk(db):
    pp = PlatePlastica.objects.create(
        description='Deep-Well Square',
        rows=8,
        columns=12,
    )
    return pp


@pytest.fixture()
def platetype_pairs(platecontext, plateplastica_pairs):
    pt = PlateType.objects.create(
        friendly='Primer pairs',
        context=platecontext,
        plastic=plateplastica_pairs,
    )
    return pt


@pytest.fixture()
def platetype_stk(platecontext, plateplastica_stk):
    pt = PlateType.objects.create(
        friendly='Primers STK',
        context=platecontext,
        plastic=plateplastica_stk,
    )
    return pt


@pytest.fixture()
def plate_united(platetype_pairs):
    p = Plate.objects.create(
        type=platetype_pairs,
        name='United_hg19_Tails_plt1',
    )
    return p


@pytest.fixture()
def plate_fwd(platetype_stk):
    p = Plate.objects.create(
        type=platetype_stk,
        name='hg19_Tails_plt1_Fw',
    )
    return p


@pytest.fixture()
def plate_rev(platetype_stk):
    p = Plate.objects.create(
        type=platetype_stk,
        name='hg19_Tails_plt1_Rev',
    )
    return p




@pytest.mark.django_db
def test_plate_united(plate_united):
    assert plate_united.name == "United_hg19_Tails_plt1"
    assert plate_united.type.plastic.columns == 12