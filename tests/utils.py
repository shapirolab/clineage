import pytest

def to_fixture(f, *args, **kwargs):
    @pytest.fixture()
    def fix(request):
        return f(request, *args, **kwargs)
    return fix
