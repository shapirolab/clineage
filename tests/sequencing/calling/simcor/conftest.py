import pytest
from sequencing.calling.models import CallingScheme


@pytest.fixture()
def simcor():
    cs = CallingScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm'
    )
    # So our objects don't have "special" objects in fields
    cs = CallingScheme.objects.get(pk=cs.pk)
    return cs
