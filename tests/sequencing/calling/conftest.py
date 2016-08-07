import pytest

from sequencing.calling.models import CallingScheme


@pytest.fixture()
def highest_peak():
    cs = CallingScheme.objects.create(
        name='highest_peak',
        description='Naive calling algorithm'
    )
    # So our objects don't have "special" objects in fields
    cs = CallingScheme.objects.get(pk=cs.pk)
    return cs