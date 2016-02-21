import pytest
import datetime

from misc.models import Taxa
from sampling.models import Individual, SampleComposition, Cell
from primers.parts.models import DNABarcode1, DNABarcode2, IlluminaReadingAdaptor1, IlluminaReadingAdaptor2
from lib_prep.workflows.models import BarcodePair, MagicalPCR1BarcodedContent, MagicalPCR1Library, AmplifiedContent
from lib_prep.multiplexes.models import PCR1MultiplexCollection, Panel
from sequencing.runs.models import MachineType, Machine, NGSKit, NGSRun
from sequencing.analysis.models import DemultiplexingScheme

from accounts.test_user import user
from misc.test_models import human_taxa
from linapp.test_models import protocoltype
from primers.parts.test_models import illuminareadingadaptor1, illuminareadingadaptor2, dnabarcode1, dnabarcode2,\
    dnabarcode1_a, dnabarcode2_a
from sampling.test_models import human_cell, human_individual, composition
from lib_prep.workflows.test_models import magicalpcr1library, pcr1multiplexcollection, magicalpcr1barcodedcontent, \
    barcodepair, barcodepair_a, magicalpcr1barcodedcontent_a,amplifiedcontent, cellcontentprotocol
from lib_prep.multiplexes.test_models import panel

@pytest.fixture()
def machinetype(db):
    mt = MachineType.objects.create(
        company="Illumina",
        model="NextSeq",
    )
    return mt


@pytest.fixture()
def machine(machinetype):
    m = Machine.objects.create(
        machineid="1",
        type=machinetype,
    )
    return m

@pytest.fixture()
def ngskit(illuminareadingadaptor1, illuminareadingadaptor2):
    nk = NGSKit.objects.create(
        reading_adaptor1=illuminareadingadaptor1,
        reading_adaptor2=illuminareadingadaptor2,
        read_length=151,
    )
    return nk


@pytest.fixture()
def ngsrun(machine, ngskit, user, magicalpcr1library, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a):
    n = NGSRun.objects.create(
        name="TestRun",
        machine=machine,
        kit=ngskit,
        user=user,
        date=datetime.date.today(),
    )
    n.libraries = [magicalpcr1library]
    return n


@pytest.mark.django_db
def test_get_samplesheet(ngsrun):
    ds = DemultiplexingScheme.objects.create(name="TestScheme", description="1")
    # TODO: sort out adaptors
    assert ngsrun.generate_sample_sheets(ds).replace('\r', '') == \
"""[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Experiment Name,TestRun,,,,,,,,
Date,{},,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,lib1,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
151,,,,,,,,,
151,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,1,,,D710,TCCGCGAA,D508,GTCAGTAC,,
2,2,,,D718,TGGGAGCC,D502,GCCTCTAT,,
""".format(datetime.date.today()).replace('\r', '')


# ------------------------------------------------------------------------------------------------

@pytest.fixture()
def full_ngsrun(db, django_user_model):
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
    b2 = DNABarcode2.objects.create(name='D508', _sequence='GTCAGTAC')
    bp12 = BarcodePair.objects.create(left=b1, right=b2)
    b3 = DNABarcode1.objects.create(name='D718', _sequence='TGGGAGCC')
    b4 = DNABarcode2.objects.create(name='D502', _sequence='GCCTCTAT')
    bp34 = BarcodePair.objects.create(left=b3, right=b4)
    b5 = DNABarcode1.objects.create(name='X720', _sequence='TGGCAGCC')
    b6 = DNABarcode2.objects.create(name='X501', _sequence='GCCTCGAT')
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
        _sequence="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    )
    ira2 = IlluminaReadingAdaptor2.objects.create(
        _sequence="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
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
    n.libraries = [mpl1, mpl2]
    return n


@pytest.mark.django_db
def test_get_samplesheet_full(full_ngsrun):
    ds = DemultiplexingScheme.objects.create(name="TestScheme", description="1")
    # TODO: sort out adaptors
    assert full_ngsrun.generate_sample_sheets(ds).replace('\r', '') == \
"""[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Experiment Name,TestRun,,,,,,,,
Date,{},,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,lib1_lib2,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
151,,,,,,,,,
151,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG,,,,,,,,
AdapterRead2,TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,1,,,D710,TCCGCGAA,D508,GTACTGAC,,
2,2,,,D718,TGGGAGCC,D502,ATAGAGGC,,
3,3,,,X720,TGGCAGCC,X501,ATCGAGGC,,
""".format(datetime.date.today()).replace('\r', '')
