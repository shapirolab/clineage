
import pytest
import datetime

from linapp.models import UserReport
from misc.models import Taxa
from sampling.models import Individual, SampleComposition, Cell
from primers.parts.models import DNABarcode1, DNABarcode2, IlluminaReadingAdaptor1, IlluminaReadingAdaptor2
from lib_prep.workflows.models import BarcodePair, MagicalPCR1BarcodedContent, MagicalPCR1Library, AmplifiedContent
from lib_prep.multiplexes.models import PCR1MultiplexCollection, Panel
from sequencing.runs.models import MachineType, Machine, NGSKit, NGSRun
from sequencing.analysis.models import DemultiplexingScheme

@pytest.fixture()
def ngsrun(db,django_user_model):
    User = django_user_model
    t = Taxa.objects.create(
        name="test_abcde",
        taxonomy_id=1,
        rank=2,
        friendly_name="test",
    )
    i = Individual.objects.create(
        taxa=t,
        sex="M",
        name="Yossi",
    )
    sc = SampleComposition.objects.create(
        name="SC"
    )
    c1 = Cell.objects.create(
        individual=i,
        name="cell1",
        composition=sc,
    )
    c2 = Cell.objects.create(
        individual=i,
        name="cell2",
        composition=sc,
    )
    b1 = DNABarcode1.objects.create(name='D710', _sequence='TCCGCGAA')
    b2 = DNABarcode2.objects.create(name='D508', _sequence='GTACTGAC')
    bp12 = BarcodePair.objects.create(left=b1, right=b2)
    b3 = DNABarcode1.objects.create(name='D718', _sequence='TGGGAGCC')
    b4 = DNABarcode2.objects.create(name='D502', _sequence='ATAGAGGC')
    bp34 = BarcodePair.objects.create(left=b3, right=b4)
    b5 = DNABarcode1.objects.create(name='X720', _sequence='TGGCAGCC')
    b6 = DNABarcode2.objects.create(name='X501', _sequence='ATCGAGGC')
    bp56 = BarcodePair.objects.create(left=b5, right=b6)
    panel = Panel.objects.create(name='mock_panel')
    pmc = PCR1MultiplexCollection.objects.create(panel=panel)
    mpl1 = MagicalPCR1Library.objects.create(name="lib1",mpx_collection=pmc)
    mpl2 = MagicalPCR1Library.objects.create(name="lib2",mpx_collection=pmc)
    ac1 = AmplifiedContent.objects.create(cell=c1,comment="1")
    ac2 = AmplifiedContent.objects.create(cell=c2,comment="2")
    mpbc1 = MagicalPCR1BarcodedContent.objects.create(barcodes=bp12, content=ac1, library=mpl1)
    mpbc2 = MagicalPCR1BarcodedContent.objects.create(barcodes=bp34, content=ac2, library=mpl1)
    mpbc3 = MagicalPCR1BarcodedContent.objects.create(barcodes=bp56, content=ac2, library=mpl2)
    mt = MachineType.objects.create(
        company="Illumina",
        model="NextSeq",
    )
    m = Machine.objects.create(
        machineid="1",
        type=mt,
    )
    ira1 = IlluminaReadingAdaptor1.objects.create(
        _sequence="AACC"
    )
    ira2 = IlluminaReadingAdaptor2.objects.create(
        _sequence="AATT"
    )
    nk = NGSKit.objects.create(
        reading_adaptor1=ira1,
        reading_adaptor2=ira2,
        read_length=151,
    )
    u = User.objects.create(
        username="Yossi",
    )
    n = NGSRun.objects.create(
        name="TestRun",
        machine=m,
        kit=nk,
        user=u,
        date=datetime.date.today(),
        #bcl_directory=,
    )
    n.libraries = [mpl1,mpl2]
    return n

@pytest.mark.django_db
def test_get_samplesheet(ngsrun):
    ds = DemultiplexingScheme.objects.create(name="TestScheme",description="1")
    assert ngsrun.generate_sample_sheets(ds) == \
"""
"""
