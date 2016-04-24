
import pytest
import tempfile
import os
import shutil

from django.conf import settings

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", help="run slow tests")


@pytest.yield_fixture(scope="session")
def temp_storage():
    t = tempfile.mkdtemp()
    settings.DATA_STORE = t
    yield
    files = os.listdir(t)
    if files:
        shutil.rmtree(t)
        raise Exception("Extra files: {}".format(files))
    else:
        os.rmdir(t)
