from sequencing.analysis.models_common import SNPHistogramGenotype, MicrosatelliteHistogramGenotype
from sequencing.runs.models import DemultiplexingScheme
from linapp.models import Protocol, ProtocolType
from lib_prep.workflows.models import CellContentProtocol


MicrosatelliteHistogramGenotype.objects.get_or_create(microsatellite=None,
        defaults=dict(repeat_number=1))
SNPHistogramGenotype.objects.get_or_create(snp=None, defaults=dict(base=""))


demultiplexingscheme, created = DemultiplexingScheme.objects.get_or_create(
    name='test demux scheme',
    description='wrovhnwpovnwecpqkewmc',
)


protocoltype, created = ProtocolType.objects.get_or_create(
    name='test protocol type'
)


cellcontentprotocol, created = CellContentProtocol.objects.get_or_create(
    initials='CCP',
    name='test cell content protocol',
    abstract='protocol abstract......\r\n.......more description',
    fulldescription='full long description......\r\n.......more description....\r\n.......',
    type=protocoltype,
)


