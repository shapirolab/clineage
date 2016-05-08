import pytest


@pytest.fixture()
def user(db,django_user_model):
    User = django_user_model
    u = User.objects.create(
        username="Shlomo",
    )
    # So our objects don't have "special" objects in fields
    u = User.objects.get(pk=u.pk)
    return u
