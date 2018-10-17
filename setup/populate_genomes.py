
from genomes.models import Chromosome, DNASlice, Assembly
from utils.ms_utils import ms_type_permutations
from targeted_enrichment.planning.models import Microsatellite
from misc.models import Taxa


human_taxa, created = Taxa.objects.get_or_create(
    name="Homo sapiens",
    taxonomy_id=9606,
    defaults=dict(
        rank='species',
        friendly_name="Human",
    )
)

hg19_assembly, created = Assembly.objects.get_or_create(
    taxa=human_taxa,
    name='February 2009 Homo sapiens (GRCh37)',
    friendly_name='hg19',
)


chrX, created = Chromosome.objects.get_or_create(
    assembly=hg19_assembly,
    name='X',
    sequence_length=155270560,
    cyclic=False,
)


chr1, created = Chromosome.objects.get_or_create(
    assembly=hg19_assembly,
    name='1',
    sequence_length=249250621,
    cyclic=False,
)