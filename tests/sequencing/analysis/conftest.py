
import pytest
import django

from distributed.utils_test import cluster
from distributed.client import Client
from tests.sequencing.conftest import *


@pytest.yield_fixture(scope="session")
def executor():
    with cluster(4) as (d, workers):
        with Client(d['address']) as c:
            def load_django():
                from django.conf import settings
                if not settings.configured:
                    databases = settings.DATABASES
                    databases['default']['NAME'] = 'test_'+databases['default']['NAME']
                    django.setup()
                django.setup()
                return settings.configured, settings.DATABASES
            c.run(load_django)
            yield c
            c.shutdown()
