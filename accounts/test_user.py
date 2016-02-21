import pytest

@pytest.fixture()
def user(db,django_user_model):
    User = django_user_model
    u = User.objects.create(
        username="Shlomo",
    )
    return u


@pytest.mark.django_db
def test_user(user):
    assert user.username == "Shlomo"
