import pytest

from wet_storage.models import PlateContext, PlatePlastica, PlateType, Plate, SampleLocation


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
