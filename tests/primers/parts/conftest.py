import pytest

from primers.parts.models import DNABarcode1, DNABarcode2, \
    IlluminaReadingAdaptor1, IlluminaReadingAdaptor2, \
    IlluminaReadingAdaptor1Cuts, IlluminaReadingAdaptor2Cuts


@pytest.fixture()
def dnabarcode1(transactional_db):
    b1 = DNABarcode1.objects.create(
        name='D508',
        _sequence='GTACTGAC'
    )
    # So our objects don't have "special" objects in fields
    b1 = DNABarcode1.objects.get(pk=b1.pk)
    return b1


@pytest.fixture()
def dnabarcode2(transactional_db):
    b2 = DNABarcode2.objects.create(
        name='D710',
        _sequence='TCCGCGAA'
    )
    # So our objects don't have "special" objects in fields
    b2 = DNABarcode2.objects.get(pk=b2.pk)
    return b2


@pytest.fixture()
def dnabarcode1_a(transactional_db):
    b3 = DNABarcode1.objects.create(
        name='D502',
        _sequence='ATAGAGGC'
    )
    # So our objects don't have "special" objects in fields
    b3 = DNABarcode1.objects.get(pk=b3.pk)
    return b3


@pytest.fixture()
def dnabarcode2_a(transactional_db):
    b4 = DNABarcode2.objects.create(
        name='D718',
        _sequence='TGGGAGCC'
    )
    # So our objects don't have "special" objects in fields
    b4 = DNABarcode2.objects.get(pk=b4.pk)
    return b4


@pytest.fixture()
def illuminareadingadaptor1(transactional_db):
    ira1 = IlluminaReadingAdaptor1.objects.create(
        name='Illumina Standard Reading Adaptor1',
        _sequence='ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    )
    # So our objects don't have "special" objects in fields
    ira1 = IlluminaReadingAdaptor1.objects.get(pk=ira1.pk)
    return ira1


@pytest.fixture()
def illuminareadingadaptor2(transactional_db):
    ira2 = IlluminaReadingAdaptor2.objects.create(
        name='Illumina Standard Reading Adaptor2',
        _sequence='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    )
    # So our objects don't have "special" objects in fields
    ira2 = IlluminaReadingAdaptor2.objects.get(pk=ira2.pk)
    return ira2


@pytest.fixture()
def illuminareadingadaptor1cuts(illuminareadingadaptor1):
    irac1 = IlluminaReadingAdaptor1Cuts.objects.create(
        ira=illuminareadingadaptor1,
        overlap_start=11,
        overlap_end=33,
    )
    # So our objects don't have "special" objects in fields
    irac1 = IlluminaReadingAdaptor1Cuts.objects.get(pk=irac1.pk)
    return irac1


@pytest.fixture()
def illuminareadingadaptor2cuts(illuminareadingadaptor2):
    irac2 = IlluminaReadingAdaptor2Cuts.objects.create(
        ira=illuminareadingadaptor2,
        overlap_start=12,
        overlap_end=34,
    )
    # So our objects don't have "special" objects in fields
    irac2 = IlluminaReadingAdaptor2Cuts.objects.get(pk=irac2.pk)
    return irac2

