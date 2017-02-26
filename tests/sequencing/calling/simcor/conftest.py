import os
import pytest
import pickle
import decimal
from sequencing.calling.models import CallingScheme
from misc.utils import get_unique_path
from sequencing.calling.models import SimulationsByCycles, FullMonoSimCorScheme, FullBiSimCorScheme, \
    ProportionalSimCorScheme, BoundProportionalSimCorScheme
from tests.sequencing.calling.conftest import *

from sequencing.calling.simcor.order.calibration.models.mutation_markov import MutationMarkov


@pytest.fixture()
def ac_model_simulations():
    # ups = [numpy.poly1d([0.00026892, -0.00205025])]
    # dws = [numpy.poly1d([0.00191615, -0.01174076]),
    #        numpy.poly1d([0.00027444, -0.00220836]),
    #        numpy.poly1d([0.0001768, -0.00199328])]
    x = [0.00026892, -0.00205025] + \
        [0.00191615, -0.01174076] + \
        [0.00027444, -0.00220836] + \
        [0.0001768, -0.00199328]

    model_params = MutationMarkov()
    model = model_params.get_for_x(x)
    sim_hists = dict()
    for cycle in range(60):
        model_cycle = model.get_for_cycles(cycle)
        for ms_len in range(3, 30):
            model_hist = model_cycle.get_hist_for_length(ms_len)
            sim_hists.setdefault(ms_len, dict())[cycle] = model_hist + ms_len
    return sim_hists


@pytest.fixture()
def simcor(ac_model_simulations):
    pickle_path = get_unique_path("pickle")
    with open(pickle_path, 'wb') as f:
        pickle.dump(ac_model_simulations, f)

    cs = SimulationsByCycles.objects.create(
        name='simulations',
        description='good old ac_mat_1a_sim_hists_15o_fresh_upto150.pickle',
        sim_hists=pickle_path,
        min_cycles=0,
        max_cycles=60,
    )
    # So our objects don't have "special" objects in fields
    cs = SimulationsByCycles.objects.get(pk=cs.pk)
    yield cs
    os.unlink(pickle_path)


@pytest.fixture()
def minimalsimcormonoschema(simcor):
    cs = FullMonoSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        min_ms_len=15,
        max_ms_len=17,
        min_cycles=20,
        max_cycles=24,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = FullMonoSimCorScheme.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def minimalsimcorbischema(simcor):
    cs = FullBiSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        min_ms_len=15,
        max_ms_len=17,
        min_cycles=20,
        max_cycles=24,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = FullBiSimCorScheme.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def minimalsimcorbipropschema(simcor):
    cs = BoundProportionalSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        lower_prop_bound=decimal.Decimal(0.4),
        upper_prop_bound=decimal.Decimal(0.6),
        min_ms_len=15,
        max_ms_len=17,
        min_cycles=20,
        max_cycles=21,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = BoundProportionalSimCorScheme.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def simcormonoschema(simcor):
    cs = FullMonoSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        min_ms_len=3,
        max_ms_len=30,
        min_cycles=20,
        max_cycles=60,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = FullMonoSimCorScheme.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def simcorbinschema(simcor):
    cs = FullBiSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        min_ms_len=3,
        max_ms_len=30,
        min_cycles=20,
        max_cycles=60,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = FullBiSimCorScheme.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def simcorbipropschema(simcor):
    cs = BoundProportionalSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        lower_prop_bound=decimal.Decimal(0.0),
        upper_prop_bound=decimal.Decimal(1.0),
        min_ms_len=3,
        max_ms_len=30,
        min_cycles=20,
        max_cycles=60,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = BoundProportionalSimCorScheme.objects.get(pk=cs.pk)
    return cs