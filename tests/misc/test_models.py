import pytest


@pytest.mark.django_db
def test_taxa(human_taxa):
    assert human_taxa.friendly_name == "Human"
