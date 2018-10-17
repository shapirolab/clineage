from primers.parts.models import DNABarcode1, DNABarcode2
from lib_prep.workflows.models import BarcodePair

dnabarcode1, created = DNABarcode1.objects.get_or_create(
    name='D508',
    _sequence='GTACTGAC'
)

dnabarcode2, created = DNABarcode2.objects.get_or_create(
    name='D710',
    _sequence='TCCGCGAA'
)

barcodepair, created = BarcodePair.objects.get_or_create(
    left=dnabarcode1,
    right=dnabarcode2
)