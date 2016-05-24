import pytest
import datetime

from lib_prep.workflows.models import BarcodePair, CellContentProtocol, AmplifiedContent, MagicalPCR1Library, MagicalPCR1BarcodedContent

from tests.primers.parts.conftest import *
from tests.sampling.conftest import *
from tests.accounts.conftest import *
from tests.linapp.conftest import *
from tests.lib_prep.multiplexes.conftest import *


@pytest.fixture()
def barcodepair(dnabarcode1, dnabarcode2):
    bp = BarcodePair.objects.create(
        left=dnabarcode1,
        right=dnabarcode2
    )
    # So our objects don't have "special" objects in fields
    bp = BarcodePair.objects.get(pk=bp.pk)
    return bp

@pytest.fixture()
def barcodepair_a(dnabarcode1_a, dnabarcode2_a):
    bp = BarcodePair.objects.create(
        left=dnabarcode1_a,
        right=dnabarcode2_a
    )
    # So our objects don't have "special" objects in fields
    bp = BarcodePair.objects.get(pk=bp.pk)
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
    # So our objects don't have "special" objects in fields
    ccp = CellContentProtocol.objects.get(pk=ccp.pk)
    return ccp


@pytest.fixture()
def amplifiedcontent(human_cell_with_se, cellcontentprotocol):
    ac = AmplifiedContent.objects.create(
        cell=human_cell_with_se,
        name='human amplified content',
        protocol=cellcontentprotocol,
        comment='some comment about this cell',
    )
    # So our objects don't have "special" objects in fields
    ac = AmplifiedContent.objects.get(pk=ac.pk)
    return ac


@pytest.fixture()
def magicalpcr1library(pcr1panel):
    mpl = MagicalPCR1Library.objects.create(
        id=1,
        name="lib1",
        panel=pcr1panel,
    )
    # So our objects don't have "special" objects in fields
    mpl = MagicalPCR1Library.objects.get(pk=mpl.pk)
    return mpl


@pytest.fixture()
def magicalpcr1barcodedcontent(barcodepair, amplifiedcontent, magicalpcr1library):
    mpbc = MagicalPCR1BarcodedContent.objects.create(
        id=1,
        barcodes=barcodepair,
        content=amplifiedcontent,
        library=magicalpcr1library,
    )
    # So our objects don't have "special" objects in fields
    mpbc = MagicalPCR1BarcodedContent.objects.get(pk=mpbc.pk)
    return mpbc


@pytest.fixture()
def magicalpcr1barcodedcontent_a(barcodepair_a, amplifiedcontent, magicalpcr1library):
    mpbc = MagicalPCR1BarcodedContent.objects.create(
        id=2,
        barcodes=barcodepair_a,
        content=amplifiedcontent,
        library=magicalpcr1library,
    )
    # So our objects don't have "special" objects in fields
    mpbc = MagicalPCR1BarcodedContent.objects.get(pk=mpbc.pk)
    return mpbc


@pytest.fixture()
def require_magicals(magicalpcr1library, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a):
    pass

