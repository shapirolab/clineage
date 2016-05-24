import pytest


@pytest.mark.django_db
def test_dnabarcode1(dnabarcode1):
    assert dnabarcode1.ref_sequence.seq == b"GTACTGAC"


@pytest.mark.django_db
def test_dnabarcode2(dnabarcode2):
    assert dnabarcode2.ref_sequence.seq == b"TCCGCGAA"


@pytest.mark.django_db
def test_illuminareadingadaptor1(illuminareadingadaptor1):
    assert illuminareadingadaptor1.ref_sequence.seq == b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"


@pytest.mark.django_db
def test_illuminareadingadaptor1(illuminareadingadaptor2):
    assert illuminareadingadaptor2.ref_sequence.seq == b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


@pytest.mark.django_db
def test_illuminareadingadaptor1fortail(illuminareadingadaptor1fortail):
    assert illuminareadingadaptor1fortail.sequence.seq == b"CTACACGACGCTCTTCCGATCT"


@pytest.mark.django_db
def test_illuminareadingadaptor2fortail(illuminareadingadaptor2fortail):
    assert illuminareadingadaptor2fortail.sequence.seq == b"CAGACGTGTGCTCTTCCGATCT"


@pytest.mark.xfail
@pytest.mark.django_db
def test_illuminareadingadaptor1forhead(illuminareadingadaptor1forhead):
    assert illuminareadingadaptor1forhead.sequence.seq == b""


@pytest.mark.xfail
@pytest.mark.django_db
def test_illuminareadingadaptor2forhead(illuminareadingadaptor2forhead):
    assert illuminareadingadaptor2forhead.sequence.seq == b""
