import pytest
import datetime

from sequencing.runs.models import MachineType, Machine, NGSKit, NGSRun, DemultiplexingScheme, Demultiplexing

from tests.primers.parts.conftest import *
from tests.lib_prep.workflows.conftest import *

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


@pytest.fixture()
def demultiplexingscheme(db):
    ds = DemultiplexingScheme.objects.create(
        name='test demux scheme',
        description='wrovhnwpovnwecpqkewmc',
    )
    return ds


@pytest.fixture()
def demultiplexing(demultiplexingscheme, ngsrun):
    dm = Demultiplexing.objects.create(
        ngs_run=ngsrun,
        demux_scheme=demultiplexingscheme
    )
    return dm
