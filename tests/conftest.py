
import pytest
import tempfile
import os

from django.conf import settings

@pytest.yield_fixture(scope="session")
def temp_storage():
    t = tempfile.mkdtemp()
    settings.DATA_STORE = t
    yield
    os.rmdir(t)
