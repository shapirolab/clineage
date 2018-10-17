from primers.parts.models import IlluminaReadingAdaptor1, IlluminaReadingAdaptor2, \
    IlluminaReadingAdaptor1ForTail, IlluminaReadingAdaptor2ForTail
from sequencing.runs.models import MachineType, Machine, NGSKit


ira1, created = IlluminaReadingAdaptor1.objects.get_or_create(
    name='Illumina Standard Reading Adaptor1',
    _sequence='ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
)


irac1, created = IlluminaReadingAdaptor1ForTail.objects.get_or_create(
    ira=ira1,
    tail_length=22,
)


ira2, created = IlluminaReadingAdaptor2.objects.get_or_create(
    name='Illumina Standard Reading Adaptor2',
    _sequence='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
)


irac2, created = IlluminaReadingAdaptor2ForTail.objects.get_or_create(
    ira=ira2,
    tail_length=22,
)


mt, created = MachineType.objects.get_or_create(
    company="Illumina",
    model="NextSeq",
    rev_left_bc=False
)


machine, created = Machine.objects.get_or_create(
    machineid="1",
    type=mt,
)


ngskit, created = NGSKit.objects.get_or_create(
    name="kit",
    reading_adaptor1=ira1,
    reading_adaptor2=ira2,
    read_length=151,
)