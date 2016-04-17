import pytest


@pytest.mark.django_db
def test_ugs(ugs_28734_right):
    assert "{}".format(ugs_28734_right) == "X:54384807-54384829(-)"


@pytest.mark.django_db
def test_te(te_28727):
    assert "{}".format(te_28727) == "X:81316094-81316116(+), X:81316243-81316265(-)"
