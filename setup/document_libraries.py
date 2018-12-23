from targeted_enrichment.amplicons.models import AmpliconCollection
from sequencing.analysis.full_msv.full_msv import get_full_ms_variations
from setup import machine, ngskit, demultiplexingscheme, amplifiedcontent, barcodepair
from sequencing.runs.models import NGSRun, Demultiplexing
from sequencing.analysis.models import SampleReads
from lib_prep.multiplexes.models import PCR1Multiplex, PCR1Panel
from lib_prep.workflows.models import MagicalPCR1Library, MagicalPCR1BarcodedContent
import datetime
from Bio import SeqIO


def document_library(amplicons, ters, fastq_r1, fastq_r2):
    ac = AmpliconCollection.custom_get_or_create(amplicons)
    fmsv = get_full_ms_variations(ac, padding=50, mss_version=0)
    pcr1m, created = PCR1Multiplex.objects.get_or_create(
        name='test PCR1Multiplex'
    )
    pcr1m.ters = ters

    pcr1panel, created = PCR1Panel.objects.get_or_create(
        name='test PCR1Panel',
        amplicon_collection=ac,
    )
    pcr1panel.mpxs = [pcr1m]

    mpl, created = MagicalPCR1Library.objects.get_or_create(
        id=1,
        name="lib1",
        panel=pcr1panel,
    )

    ngsrun, created = NGSRun.objects.get_or_create(
        name="TestRun",
        machine=machine,
        kit=ngskit,
        date=datetime.date.today(),
        bcl_directory=''
    )
    ngsrun.libraries = [mpl]

    demultiplexing, created = Demultiplexing.objects.get_or_create(
        ngs_run=ngsrun,
        demux_scheme=demultiplexingscheme
    )

    bc, created = MagicalPCR1BarcodedContent.objects.get_or_create(
        id=1,
        barcodes=barcodepair,
        content=amplifiedcontent,
        library=mpl,
    )
    with open(fastq_r1) as f:
        seq = SeqIO.parse(f, format="fastq")
        num_reads = sum(1 for _ in seq)

    sr = SampleReads.objects.create(
        demux=demultiplexing,
        barcoded_content=bc,
        library=mpl,
        fastq1=fastq_r1,
        fastq2=fastq_r2,
        num_reads=num_reads,
    )

    return sr