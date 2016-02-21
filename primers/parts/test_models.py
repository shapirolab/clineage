import pytest
from models import DNABarcode1, DNABarcode2


@pytest.fixture()
def dnabarcode1(db):
    b1 = DNABarcode1.objects.create(
        name='D710',
        _sequence='TCCGCGAA'
    )
    return b1


@pytest.fixture()
def dnabarcode2(db):
    b2 = DNABarcode2.objects.create(
        name='D508',
        _sequence='GTACTGAC'
    )
    return b2


@pytest.fixture()
def dnabarcode1_a(db):
    b3 = DNABarcode1.objects.create(
        name='D718',
        _sequence='TGGGAGCC'
    )
    return b3


@pytest.fixture()
def dnabarcode2_a(db):
    b4 = DNABarcode2.objects.create(
        name='D502',
        _sequence='ATAGAGGC'
    )
    return b4


@pytest.mark.django_db
def test_dnabarcode1(dnabarcode1):
    assert dnabarcode1._sequence == "TCCGCGAA"


@pytest.mark.django_db
def test_dnabarcode2(dnabarcode2):
    assert dnabarcode2._sequence == "GTACTGAC"