from sampling.models import Individual, Cell, SampleComposition
from setup.populate_genomes import human_taxa
from setup.populate_misc import cellcontentprotocol
from django.contrib.auth.models import User
from django.contrib.auth.hashers import make_password
from lib_prep.workflows.models import AmplifiedContent
from lib_prep.workflows.models import MagicalPCR1BarcodedContent

sc, created = SampleComposition.objects.get_or_create(
    name="Single Cell"
)


user, created = User.objects.get_or_create(
    username="Shlomo",
    defaults=dict(
        password=make_password("abcd"),
    )
)


human_individual, created = Individual.objects.get_or_create(
    taxa=human_taxa,
    sex="M",
    name="Yossi",
    partner=user,
)


cell, created = Cell.objects.get_or_create(
    individual=human_individual,
    name='human_cell_no_se',
    composition=sc,
)


amplifiedcontent, created = AmplifiedContent.objects.get_or_create(
    cell=cell,
    name='human amplified content',
    protocol=cellcontentprotocol,
    comment='some comment about this cell',
)