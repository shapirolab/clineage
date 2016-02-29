import pytest


@pytest.mark.django_db
def test_te(te_28727):
    assert te_28727.targets.get(name='X_81316201_81316236').slice.start_pos == 81316201
