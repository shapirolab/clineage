import pytest

from misc.models import Taxa


@pytest.fixture()
def human_taxa(transactional_db):
    t = Taxa.objects.create(
        name="Homo sapiens",
        taxonomy_id=9606,
        rank='species',
        friendly_name="Human",
    )
    return t
