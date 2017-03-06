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
def hg19_chromosome1(hg19_assembly):
    c = Chromosome.objects.create(
        assembly=hg19_assembly,
        name='1',
        sequence_length=249250621,
        cyclic=False,
    )
    # So our objects don't have "special" objects in fields
    c = Chromosome.objects.get(pk=c.pk)
    return c


@pytest.fixture()
def hg19_chromosome4(hg19_assembly):
    c = Chromosome.objects.create(
        assembly=hg19_assembly,
        name='4',
        sequence_length=191154276,
        cyclic=False,
    )
    # So our objects don't have "special" objects in fields
    c = Chromosome.objects.get(pk=c.pk)
    return c


@pytest.fixture()
def hg19_chromosome10(hg19_assembly):
    c = Chromosome.objects.create(
        assembly=hg19_assembly,
        name='10',
        sequence_length=135534747,
        cyclic=False,
    )
    # So our objects don't have "special" objects in fields
    c = Chromosome.objects.get(pk=c.pk)
    return c


@pytest.fixture()
def hg19_chromosome21(hg19_assembly):
    c = Chromosome.objects.create(
        assembly=hg19_assembly,
        name='21',
        sequence_length=48129895,
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
        start_pos=39873585,
        end_pos=39873604,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_target_a(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=39873639,
        end_pos=39873678,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_target_b(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=39873679,
        end_pos=39873696,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_right(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=39873745,
        end_pos=39873766,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_adj_ms_2_amplicon(hg19_chromosome):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome,
        start_pos=39873605,
        end_pos=39873744,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_1_left(hg19_chromosome1):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome1,
        start_pos=115256468,
        end_pos=115256487,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_1_target_a(hg19_chromosome1):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome1,
        start_pos=115256524,
        end_pos=115256524,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_1_target_b(hg19_chromosome1):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome1,
        start_pos=115256529,
        end_pos=115256529,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_1_right(hg19_chromosome1):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome1,
        start_pos=115256593,
        end_pos=115256612,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_1_amplicon(hg19_chromosome1):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome1,
        start_pos=115256488,
        end_pos=115256592,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_2_left(hg19_chromosome10):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome10,
        start_pos=10341329,
        end_pos=10341354,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_2_target_a(hg19_chromosome10):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome10,
        start_pos=10341467,
        end_pos=10341467,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_2_right(hg19_chromosome10):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome10,
        start_pos=10341543,
        end_pos=10341565,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_2_amplicon(hg19_chromosome10):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome10,
        start_pos=10341355,
        end_pos=10341542,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_3_left(hg19_chromosome4):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome4,
        start_pos=106155054,
        end_pos=106155076,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_3_right(hg19_chromosome4):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome4,
        start_pos=106155241,
        end_pos=106155263,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_3_amplicon(hg19_chromosome4):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome4,
        start_pos=106155077,
        end_pos=106155240,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_4_left(hg19_chromosome21):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome21,
        start_pos=36259181,
        end_pos=36259202,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_4_right(hg19_chromosome21):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome21,
        start_pos=36259420,
        end_pos=36259442,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas


@pytest.fixture()
def slice_snp_4_amplicon(hg19_chromosome21):
    dnas = DNASlice.objects.create(
        chromosome=hg19_chromosome21,
        start_pos=36259203,
        end_pos=36259419,
    )
    # So our objects don't have "special" objects in fields
    dnas = DNASlice.objects.get(pk=dnas.pk)
    return dnas
