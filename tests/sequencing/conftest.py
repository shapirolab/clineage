import pytest
from sequencing.analysis.models_common import SNPHistogramGenotype, MicrosatelliteHistogramGenotype


@pytest.fixture()
def requires_none_genotypes(request, transactional_db):
    MicrosatelliteHistogramGenotype.objects.get_or_create(microsatellite=None,
        defaults=dict(repeat_number=1))
    SNPHistogramGenotype.objects.get_or_create(snp=None, defaults=dict(base=""))