import pytest

from linapp.models import Protocol, ProtocolType


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
