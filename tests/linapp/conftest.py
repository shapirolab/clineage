import pytest
import os

from linapp.models import Protocol, ProtocolType
from tests.sampling.conftest import *


@pytest.fixture()
def protocoltype(transactional_db):
    pt = ProtocolType.objects.create(
        name='test protocol type'
    )
    # So our objects don't have "special" objects in fields
    pt = ProtocolType.objects.get(pk=pt.pk)
    return pt


@pytest.fixture()
def protocol(protocoltype):
    p = Protocol.objects.create(
        initials='TP',
        name='test protocol',
        abstract='protocol abstract......\r\n.......more description',
        fulldescription='full long description......\r\n.......more description....\r\n.......',
        type=protocoltype,
    )
    # So our objects don't have "special" objects in fields
    p = Protocol.objects.get(pk=p.pk)
    return p

@pytest.fixture()
def cell_file_no_se(human_cell_no_se):
    cell_name=human_cell_no_se.id
    partner=human_cell_no_se
    file="{}-{}.hist".format(cell_name, partner)
    print("cell_file_no_se dir:")
    print(os.getcwd())
    with open(file, 'wb') as f:
        pass
    # yield (file)
    return file
    os.unlink(file)

@pytest.fixture()
def cell_file_with_se(human_cell_with_se):
    cell_name=human_cell_with_se.id
    partner=human_cell_with_se
    file="{}-{}.hist".format(cell_name, partner)
    print("cell_file_with_se dir:")
    print(os.getcwd())
    with open(file, 'wb') as f:
        pass
    # yield (file)
    return file
    os.unlink(file)