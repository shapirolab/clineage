import pytest
from models import BarcodePair, CellContentProtocol, AmplifiedContent, MagicalPCR1Library, MagicalPCR1BarcodedContent
from ..multiplexes.test_models import pcr1multiplexcollection
from linapp.test_models import protocoltype
from primers.parts.test_models import dnabarcode1, dnabarcode2
from sampling.test_models import human_cell, human_individual, human_taxa, composition, human_samplingevent, human_extraction, human_extractionevent, organ, tissue
from accounts.test_user import user
import datetime


@pytest.fixture()
def barcodepair(dnabarcode1, dnabarcode2):
    bp = BarcodePair.objects.create(
        left=dnabarcode1,
        right=dnabarcode2
    )
    return bp


@pytest.fixture()
def cellcontentprotocol(protocoltype):
    ccp = CellContentProtocol.objects.create(
        initials='CCP',
        name='test cell content protocol',
        abstract='protocol abstract......\r\n.......more description',
        fulldescription='full long description......\r\n.......more description....\r\n.......',
        type=protocoltype,
    )
    return ccp


@pytest.fixture()
def amplifiedcontent(human_cell, cellcontentprotocol):
    ac = AmplifiedContent.objects.create(
        cell=human_cell,
        name='human amplified content',
        protocol=cellcontentprotocol,
        comment='some comment about this cell',
    )
    return ac


@pytest.fixture()
def magicalpcr1library(pcr1multiplexcollection):
    mpl = MagicalPCR1Library.objects.create(
        name="lib1",
        mpx_collection=pcr1multiplexcollection,
    )
    return mpl


@pytest.fixture()
def magicalpcr1barcodedcontent(barcodepair, amplifiedcontent, magicalpcr1library):
    mpbc = MagicalPCR1BarcodedContent.objects.create(
        barcodes=barcodepair,
        content=amplifiedcontent,
        library=magicalpcr1library,
    )
    return mpbc


@pytest.mark.django_db
def test_amplifiedcontent(amplifiedcontent):
    assert amplifiedcontent.name == 'human amplified content'
