import decimal
import pickle
from misc.utils import get_unique_path
from sequencing.calling.simcor.order.calibration.models.mutation_markov import MutationMarkov
from sequencing.calling.models import HighestPeaksProximityRatioFilteredBiSimCorSchemeModel, SimulationsByCycles
from sequencing.calling.models import HighestPeaksMonoSimCorSchemeModel


def ac_model_simulations():
    x = [0.00026892, -0.00205025] + \
        [0.00191615, -0.01174076] + \
        [0.00027444, -0.00220836] + \
        [0.0001768, -0.00199328]

    model_params = MutationMarkov()
    model = model_params.get_for_x(x)
    sim_hists = dict()
    for cycle in range(60):
        model_cycle = model.get_for_cycles(cycle)
        for ms_len in range(3, 40):
            model_hist = model_cycle.get_hist_for_length(ms_len)
            sim_hists.setdefault(ms_len, dict())[cycle] = model_hist + ms_len
    return sim_hists


def polyg_model_simulations():
    x = [-0.00089469, 0.00538763] + \
        [0.00496469, -0.0378956] + \
        [0.00127528, -0.01155447] + \
        [0.00029054, -0.00350671]

    model_params = MutationMarkov()
    model = model_params.get_for_x(x)
    sim_hists = dict()
    for cycle in range(60):
        model_cycle = model.get_for_cycles(cycle)
        for ms_len in range(3, 20):
            model_hist = model_cycle.get_hist_for_length(ms_len)
            sim_hists.setdefault(ms_len, dict())[cycle] = model_hist + ms_len
    return sim_hists


def polya_model_simulations():
    x = [0.0012521, -0.0091638] + \
        [0.00408945, -0.03640044] + \
        [0.00067988, -0.00947355] + \
        [0.00025592, -0.00343027]

    model_params = MutationMarkov()
    model = model_params.get_for_x(x)
    sim_hists = dict()
    for cycle in range(60):
        model_cycle = model.get_for_cycles(cycle)
        for ms_len in range(3, 45):
            model_hist = model_cycle.get_hist_for_length(ms_len)
            sim_hists.setdefault(ms_len, dict())[cycle] = model_hist + ms_len
    return sim_hists


pickle_path = get_unique_path("pickle")
with open(pickle_path, 'wb') as f:
    pickle.dump(ac_model_simulations(), f)

sbc_ac, created = SimulationsByCycles.objects.get_or_create(
    name='AC markov simulations',
    description='AC simulations 0-60 PCR cycles',
    defaults=dict(
        sim_hists=pickle_path,
        min_cycles=0,
        max_cycles=60,
    )
)

ac_schema_bi, created = HighestPeaksProximityRatioFilteredBiSimCorSchemeModel.objects.get_or_create(
    name='AC markov HP-PRF strinct extended fast',
    description='Simulations correlation highest peaks calling algorithm',
    proportion_step=decimal.Decimal(0.1),
    lower_prop_bound=decimal.Decimal(0.1),
    upper_prop_bound=decimal.Decimal(0.9),
    min_ms_len=5,
    max_ms_len=30,
    min_cycles=10,
    max_cycles=55,
    range_from_point=2,
    minimal_seeds_distance=3,
    length_sensitivity=decimal.Decimal(0.25),  # Correction following exclusion function revision (removing steps)
    diff_sensetivity=decimal.Decimal(0.75),  # Correction following exclusion function revision (removing steps)
    cycle_sensetivity=decimal.Decimal(1.0),  # This is currently doing nothing
    simulations=sbc_ac
)


pickle_path = get_unique_path("pickle")
with open(pickle_path, 'wb') as f:
    pickle.dump(polyg_model_simulations(), f)

sbc_g, created = SimulationsByCycles.objects.get_or_create(
    name='G markov simulations',
    description='PolyG simulations 0-60 PCR cycles',
    defaults=dict(
        sim_hists=pickle_path,
        min_cycles=0,
        max_cycles=60,
    ),
)


pickle_path = get_unique_path("pickle")
with open(pickle_path, 'wb') as f:
    pickle.dump(polya_model_simulations(), f)

sbc_a, created = SimulationsByCycles.objects.get_or_create(
    name='A markov simulations',
    description='PolyG simulations 0-60 PCR cycles',
    defaults=dict(
        sim_hists=pickle_path,
        min_cycles=0,
        max_cycles=60,
    ),
)


a_schema_mono, created = HighestPeaksMonoSimCorSchemeModel.objects.get_or_create(
    name='PolyA markov mono',
    description='Simulations correlation highest peaks calling algorithm',
    min_ms_len=3,
    max_ms_len=30,
    min_cycles=5,
    max_cycles=60,
    range_from_point=3,
    simulations=sbc_a
)
g_schema_mono, created = HighestPeaksMonoSimCorSchemeModel.objects.get_or_create(
    name='PolyG markov mono',
    description='Simulations correlation highest peaks calling algorithm',
    min_ms_len=3,
    max_ms_len=20,
    min_cycles=5,
    max_cycles=60,
    range_from_point=3,
    simulations=sbc_g
)
ac_schema_mono, created = HighestPeaksMonoSimCorSchemeModel.objects.get_or_create(
    name='AC markov mono',
    description='Simulations correlation highest peaks calling algorithm',
    min_ms_len=5,
    max_ms_len=30,
    min_cycles=5,
    max_cycles=60,
    range_from_point=2,
    simulations=sbc_ac
)