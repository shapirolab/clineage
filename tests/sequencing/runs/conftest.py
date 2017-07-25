import pytest
import datetime
import os

from sequencing.runs.models import MachineType, Machine, NGSKit, NGSRun, DemultiplexingScheme, \
    Demultiplexing, MergedDemultiplexing

from tests.primers.parts.conftest import *
from tests.lib_prep.workflows.conftest import *

from misc.utils import get_unique_path

@pytest.fixture()
def machinetype(transactional_db):
    mt = MachineType.objects.create(
        company="Illumina",
        model="NextSeq",
        rev_left_bc=False
    )
    # So our objects don't have "special" objects in fields
    mt = MachineType.objects.get(pk=mt.pk)
    return mt


@pytest.fixture()
def machine(machinetype):
    m = Machine.objects.create(
        machineid="1",
        type=machinetype,
    )
    # So our objects don't have "special" objects in fields
    m = Machine.objects.get(pk=m.pk)
    return m

@pytest.fixture()
def ngskit(illuminareadingadaptor1, illuminareadingadaptor2):
    nk = NGSKit.objects.create(
        name="kit",
        reading_adaptor1=illuminareadingadaptor1,
        reading_adaptor2=illuminareadingadaptor2,
        read_length=151,
    )
    # So our objects don't have "special" objects in fields
    nk = NGSKit.objects.get(pk=nk.pk)
    return nk


@pytest.yield_fixture()
def ngsrun(machine, ngskit, magicalpcr1library, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a, temp_storage):
    d = get_unique_path()
    os.mkdir(d)
    n = NGSRun.objects.create(
        name="TestRun",
        machine=machine,
        kit=ngskit,
        date=datetime.date.today(),
        bcl_directory=d
    )
    n.libraries = [magicalpcr1library]
    # So our objects don't have "special" objects in fields
    n = NGSRun.objects.get(pk=n.pk)
    yield n
    os.rmdir(d)

@pytest.yield_fixture()
def ngsrun2(machine, ngskit, magicalpcr1library, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a, temp_storage):
    d = get_unique_path()
    os.mkdir(d)
    n = NGSRun.objects.create(
        name="Second TestRun",
        machine=machine,
        kit=ngskit,
        date=datetime.date.today(),
        bcl_directory=d
    )
    n.libraries = [magicalpcr1library]
    # So our objects don't have "special" objects in fields
    n = NGSRun.objects.get(pk=n.pk)
    yield n
    os.rmdir(d)

@pytest.fixture()
def demultiplexingscheme(transactional_db):
    ds = DemultiplexingScheme.objects.create(
        name='test demux scheme',
        description='wrovhnwpovnwecpqkewmc',
    )
    # So our objects don't have "special" objects in fields
    ds = DemultiplexingScheme.objects.get(pk=ds.pk)
    return ds


@pytest.fixture()
def demultiplexing(demultiplexingscheme, ngsrun):
    dm = Demultiplexing.objects.create(
        ngs_run=ngsrun,
        demux_scheme=demultiplexingscheme
    )
    # So our objects don't have "special" objects in fields
    dm = Demultiplexing.objects.get(pk=dm.pk)
    return dm


@pytest.fixture()
def mergeddemultiplexing(demultiplexingscheme, ngsrun, ngsrun2):
    dm = MergedDemultiplexing.objects.create(
        ngs_run=ngsrun,
        demux_scheme=demultiplexingscheme,
    )
    dm.ngs_runs = [ngsrun2]
    dm.save()
    # So our objects don't have "special" objects in fields
    dm = Demultiplexing.objects.get(pk=dm.pk)
    return dm
