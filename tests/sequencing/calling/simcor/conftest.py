import os
import pytest
import pickle
import decimal
from sequencing.calling.models import CallingScheme
from misc.utils import get_unique_path
from sequencing.calling.models import SimulationsByCycles, FullMonoSimCorScheme, FullBiSimCorScheme, \
    ProportionalSimCorScheme, BoundProportionalSimCorScheme, HighestPeaksProportionalBiSimCorSchemeModel, \
    ProximityRatioFilteredBoundProportionalSimCorScheme, HighestPeaksProximityRatioFilteredBiSimCorSchemeModel,\
    HighestPeaksMonoSimCorSchemeModel
from tests.sequencing.calling.conftest import *

from sequencing.calling.simcor.order.calibration.models.mutation_markov import MutationMarkov
from sequencing.calling.hist import Histogram as dHist


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


@pytest.yield_fixture()
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
    cs = ProportionalSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        proportion_step=decimal.Decimal(0.1),
        min_ms_len=15,
        max_ms_len=17,
        min_cycles=20,
        max_cycles=21,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = ProportionalSimCorScheme.objects.get(pk=cs.pk)
    return cs



@pytest.fixture()
def minimalsimcorbiboundpropschema(simcor):
    cs = BoundProportionalSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        proportion_step=decimal.Decimal(0.1),
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
        proportion_step=decimal.Decimal(0.1),
        lower_prop_bound=decimal.Decimal(0.3),
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


@pytest.fixture()
def minimal_prf_simcorbiboundpropschema(simcor):
    cs = ProximityRatioFilteredBoundProportionalSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        proportion_step=decimal.Decimal(0.1),
        lower_prop_bound=decimal.Decimal(0.4),
        upper_prop_bound=decimal.Decimal(0.6),
        min_ms_len=15,
        max_ms_len=17,
        min_cycles=20,
        max_cycles=21,
        length_sensitivity=decimal.Decimal(0.11),  # Correction following exclusion function revision (removing steps)
        diff_sensetivity=decimal.Decimal(0.75),  # Correction following exclusion function revision (removing steps)
        cycle_sensetivity=decimal.Decimal(1.0),
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = ProximityRatioFilteredBoundProportionalSimCorScheme.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def prf_simcorbipropschema(simcor):
    cs = ProximityRatioFilteredBoundProportionalSimCorScheme.objects.create(
        name='simcor',
        description='Simulations correlation calling algorithm',
        proportion_step=decimal.Decimal(0.1),
        lower_prop_bound=decimal.Decimal(0.3),
        upper_prop_bound=decimal.Decimal(1.0),
        min_ms_len=3,
        max_ms_len=30,
        min_cycles=20,
        max_cycles=60,
        length_sensitivity=decimal.Decimal(0.11),  # Correction following exclusion function revision (removing steps)
        diff_sensetivity=decimal.Decimal(0.75),  # Correction following exclusion function revision (removing steps)
        cycle_sensetivity=decimal.Decimal(1.0),
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = ProximityRatioFilteredBoundProportionalSimCorScheme.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def simcorbiprophighpeakschema(simcor):
    cs = HighestPeaksProportionalBiSimCorSchemeModel.objects.create(
        name='simcor',
        description='Simulations correlation highest peaks calling algorithm',
        proportion_step=decimal.Decimal(0.1),
        lower_prop_bound=decimal.Decimal(0.3),
        upper_prop_bound=decimal.Decimal(1.0),
        min_ms_len=3,
        max_ms_len=30,
        min_cycles=20,
        max_cycles=60,
        range_from_point=3,
        minimal_seeds_distance=3,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = HighestPeaksProportionalBiSimCorSchemeModel.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def prf_simcorbiprophighpeakschema(simcor):
    cs = HighestPeaksProximityRatioFilteredBiSimCorSchemeModel.objects.create(
        name='simcor',
        description='Simulations correlation highest peaks calling algorithm',
        proportion_step=decimal.Decimal(0.1),
        lower_prop_bound=decimal.Decimal(0.3),
        upper_prop_bound=decimal.Decimal(1.0),
        min_ms_len=3,
        max_ms_len=30,
        min_cycles=20,
        max_cycles=60,
        range_from_point=3,
        minimal_seeds_distance=3,
        length_sensitivity=decimal.Decimal(0.21),
        diff_sensetivity=decimal.Decimal(0.65),
        cycle_sensetivity=decimal.Decimal(1.0),
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = HighestPeaksProximityRatioFilteredBiSimCorSchemeModel.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def minimal_prf_simcorbiprophighpeakschema(simcor):
    cs = HighestPeaksProximityRatioFilteredBiSimCorSchemeModel.objects.create(
        name='simcor',
        description='Simulations correlation highest peaks calling algorithm',
        proportion_step=decimal.Decimal(0.1),
        lower_prop_bound=decimal.Decimal(0.4),
        upper_prop_bound=decimal.Decimal(0.6),
        min_ms_len=15,
        max_ms_len=17,
        min_cycles=20,
        max_cycles=21,
        range_from_point=3,
        minimal_seeds_distance=3,
        length_sensitivity=decimal.Decimal(0.11),  # Correction following exclusion function revision (removing steps)
        diff_sensetivity=decimal.Decimal(0.75),  # Correction following exclusion function revision (removing steps)
        cycle_sensetivity=decimal.Decimal(1.0),
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = HighestPeaksProximityRatioFilteredBiSimCorSchemeModel.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def simcormonoprophighpeakschema(simcor):
    cs = HighestPeaksMonoSimCorSchemeModel.objects.create(
        name='simcor',
        description='Simulations correlation highest peaks calling algorithm',
        min_ms_len=3,
        max_ms_len=30,
        min_cycles=20,
        max_cycles=60,
        range_from_point=3,
        simulations=simcor
    )
    # So our objects don't have "special" objects in fields
    cs = HighestPeaksMonoSimCorSchemeModel.objects.get(pk=cs.pk)
    return cs


@pytest.fixture()
def HP_HISTOGRAMS():
    HP_HISTOGRAM_1 = (6,15),{
    6:60,
    7:59,
    8:58,
    9:57,
    10:5,
    11:4,
    12:57,
    13:58,
    14:59,
    15:60,
        }

    HP_HISTOGRAM_2 = (15,7), {
        6:8,
        7:59,
        8:58,
        9:57,
        10:5,
        11:4,
        12:57,
        13:58,
        14:59,
        15:60,
    }

    HP_HISTOGRAM_3 = (11,6), {
        6:60,
        7:59,
        8:58,
        9:57,
        10:5,
        11:150,
        12:57,
        13:58,
        14:59,
        15:59,
    }

    HP_HISTOGRAM_4 = (11,15), {
        6:59,
        7:59,
        8:58,
        9:57,
        10:5,
        11:150,
        12:57,
        13:58,
        14:59,
        15:60,
    }

    HP_HISTOGRAM_5 = (8,15), {
        6:61,
        7:59,
        8:63,
        9:57,
        10:5,
        11:4,
        12:57,
        13:58,
        14:59,
        15:60,
    }
    return HP_HISTOGRAM_1, HP_HISTOGRAM_2, HP_HISTOGRAM_3, HP_HISTOGRAM_4, HP_HISTOGRAM_5


@pytest.fixture()
def histogram_dicts(HP_HISTOGRAMS):
    HP_HISTOGRAM_1, HP_HISTOGRAM_2, HP_HISTOGRAM_3, HP_HISTOGRAM_4, HP_HISTOGRAM_5 = HP_HISTOGRAMS
    hp1 = HP_HISTOGRAM_1[0], dHist(HP_HISTOGRAM_1[1])
    hp2 = HP_HISTOGRAM_2[0], dHist(HP_HISTOGRAM_2[1])
    hp3 = HP_HISTOGRAM_3[0], dHist(HP_HISTOGRAM_3[1])
    hp4 = HP_HISTOGRAM_4[0], dHist(HP_HISTOGRAM_4[1])
    hp5 = HP_HISTOGRAM_5[0], dHist(HP_HISTOGRAM_5[1])
    hists = {
        1: {'result': hp1[0], 'histogram': hp1[1]},
        2: {'result': hp2[0], 'histogram': hp2[1]},
        3: {'result': hp3[0], 'histogram': hp3[1]},
        4: {'result': hp4[0], 'histogram': hp4[1]},
        5: {'result': hp5[0], 'histogram': hp5[1]},
    }
    return hists
