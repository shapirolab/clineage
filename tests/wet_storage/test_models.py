import pytest


@pytest.mark.django_db
def test_plate_united(plate_united):
    assert plate_united.name == "United_hg19_Tails_plt1"
    assert plate_united.type.plastic.columns == 12
