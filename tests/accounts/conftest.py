import pytest
from django.contrib.auth.hashers import make_password


@pytest.fixture()
def user(db,django_user_model):
    User = django_user_model
    u = User.objects.create(
        username="Shlomo",
        password=make_password("abcd"),
    )
    # So our objects don't have "special" objects in fields
    u = User.objects.get(pk=u.pk)
    return u

@pytest.fixture()
def loggedin_client(client):
    resp = client.post('/accounts/login/', {'username': 'Shlomo', 'password': 'abcd'})
    assert resp.status_code in [200, 302]
    return client
