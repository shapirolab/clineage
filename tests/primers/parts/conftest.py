import pytest

from primers.parts.models import DNABarcode1, DNABarcode2, \
    IlluminaReadingAdaptor1, IlluminaReadingAdaptor2, \
    IlluminaReadingAdaptor1Cuts, IlluminaReadingAdaptor2Cuts


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


@pytest.fixture()
def illuminareadingadaptor1(db):
    ira1 = IlluminaReadingAdaptor1.objects.create(
        name='Illumina Standard Reading Adaptor1',
        _sequence='ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    )
    return ira1


@pytest.fixture()
def illuminareadingadaptor2(db):
    ira2 = IlluminaReadingAdaptor2.objects.create(
        name='Illumina Standard Reading Adaptor2',
        _sequence='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    )
    return ira2


@pytest.fixture()
def illuminareadingadaptor1cuts(illuminareadingadaptor1):
    irac1 = IlluminaReadingAdaptor1Cuts.objects.create(
        ira=illuminareadingadaptor1,
        overlap_start=11,
        overlap_end=33,
    )
    return irac1


@pytest.fixture()
def illuminareadingadaptor2cuts(illuminareadingadaptor2):
    irac2 = IlluminaReadingAdaptor2Cuts.objects.create(
        ira=illuminareadingadaptor2,
        overlap_start=12,
        overlap_end=34,
    )
    return irac2

