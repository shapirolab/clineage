import pytest
from models import Taxa

@pytest.fixture()
def human_taxa(db):
    t = Taxa.objects.create(
        name="Homo sapiens",
        taxonomy_id=9606,
        rank='species',
        friendly_name="Human",
    )
    return t

@pytest.mark.django_db
def test_taxa(human_taxa):
    assert human_taxa.friendly_name == "Human"