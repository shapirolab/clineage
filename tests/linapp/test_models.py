import pytest


@pytest.mark.django_db
def test_protocol(protocol):
    assert protocol.initials == "TP"
