import pytest

from genomes.models import Assembly, Chromosome, DNASlice

from tests.misc.conftest import human_taxa


@pytest.fixture()
def hg19_assembly(human_taxa):
    a = Assembly.objects.create(
        taxa=human_taxa,
        name='February 2009 Homo sapiens (GRCh37)',
        friendly_name='hg19',
    )
    # So our objects don't have "special" objects in fields
    a = Assembly.objects.get(pk=a.pk)
    return a


@pytest.fixture()
def hg19_chromosome(hg19_assembly):
    c = Chromosome.objects.create(
        assembly=hg19_assembly,
        name='X',
        sequence_length=155270560,
        cyclic=False,
    )
    # So our objects don't have "special" objects in fields
    c = Chromosome.objects.get(pk=c.pk)
    return c


@pytest.fixture()
def slice_28727_left(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=81316094,
        end_pos=81316116,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28727_target_a(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=81316131,
        end_pos=81316199,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28727_target_b(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=81316201,
        end_pos=81316236,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28727_right(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=81316243,
        end_pos=81316265,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28727_amplicon(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=81316117,
        end_pos=81316242,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28727_overlaps_some(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=81316110,
        end_pos=81316142,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28734_left(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=54384674,
        end_pos=54384696,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28734_target_a(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=54384788,
        end_pos=54384805,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28734_right(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=54384807,
        end_pos=54384829,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_28734_amplicon(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=54384697,
        end_pos=54384806,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_1_left(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=74123138,
        end_pos=74123160,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_1_target_a(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=74123161,
        end_pos=74123220,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_1_right(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=74123282,
        end_pos=74123304,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_1_amplicon(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=74123161,
        end_pos=74123281,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_left(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=35580142,
        end_pos=35580164,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_target_a(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=35580220,
        end_pos=35580271,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_target_b(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=35580272,
        end_pos=35580315,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_right(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=35580342,
        end_pos=35580364,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_amplicon(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=35580165,
        end_pos=35580341,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas
