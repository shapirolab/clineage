import pytest


@pytest.mark.django_db
def test_user(user):
    assert user.username == "Shlomo"
