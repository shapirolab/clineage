import pytest
from models import Protocol, ProtocolType


@pytest.fixture()
def protocoltype(db):
    pt = ProtocolType.objects.create(
        name='test protocol type'
    )
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
    return p


@pytest.mark.django_db
def test_protocol(protocol):
    assert protocol.initials == "TP"
