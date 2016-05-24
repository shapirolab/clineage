import pytest


@pytest.mark.django_db
def test_individual(human_individual):
    assert human_individual.name == "Yossi"


@pytest.mark.django_db
def test_cell_without_sampling(human_cell_no_se):
    assert human_cell_no_se.name == "human_cell_no_se"
    assert human_cell_no_se.sampling is None


@pytest.mark.django_db
def test_cell_with_sampling(human_cell_with_se):
    assert human_cell_with_se.name == "human_cell_with_se"
    assert human_cell_with_se.sampling.name == "human_samplingevent"
