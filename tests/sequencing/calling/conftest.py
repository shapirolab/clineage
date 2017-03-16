import pytest
from sequencing.analysis.full_msv.models import SampleReads, FullMSVMergedReads, \
    FullMSVAssignment, FullMSVHistogram, FullMSVariations, Histogram
from sequencing.analysis.models_common import PearOutputMixin, SNPHistogramGenotypeSet, \
    MicrosatelliteHistogramGenotypeSet, SNPHistogramGenotype, MicrosatelliteHistogramGenotype, \
    HistogramEntryReads
from misc.utils import get_unique_path

from tests.sequencing.conftest import *
from tests.sequencing.runs.conftest import *
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.amplicons.conftest import *
from tests.sequencing.calling.hists_dict import MS_HISTOGRAMS_DICT
from tests.targeted_enrichment.planning.conftest import *
from tests.flat_dict import FlatDict
from tests.sequencing.analysis.conftest import *


def touch(fp):
    with open(fp, 'wb') as f:
        pass


@pytest.fixture(scope="session")
def histograms_fd():
    return MS_HISTOGRAMS_DICT


@pytest.fixture()
def dummy_samplereads(demultiplexing, magicalpcr1barcodedcontent):
    """
    Generate dummy SampleReads objects with empty fastq files
    """
    fastq_r1 = get_unique_path("fastq")
    touch(fastq_r1)
    fastq_r2 = get_unique_path("fastq")
    touch(fastq_r2)
    sr = SampleReads.objects.create(
        demux=demultiplexing,
        barcoded_content=magicalpcr1barcodedcontent,
        library=demultiplexing.ngs_run.libraries.get(),
        num_reads=100,  # Dummy
        fastq1=fastq_r1,
        fastq2=fastq_r2,
    )
    # So our objects don't have "special" objects in fields
    sr = SampleReads.objects.get(pk=sr.pk)
    yield sr
    sr.delete()


@pytest.fixture()
def histograms_and_calling_solutions_d(dummy_samplereads, histograms_fd, requires_microsatellites, requires_none_genotypes):
    none_snp_genotype = SNPHistogramGenotype.objects.get(snp=None)
    snp_histogram_genotypes, c = SNPHistogramGenotypeSet.objects.get_or_create(
        **{fn: none_snp_genotype for fn in SNPHistogramGenotypeSet.genotype_field_names()})
    d = dict()
    for amp_id, ms_dict in histograms_fd.items():
        d[amp_id] = dict()
        for ms_id, genotypes_dict in ms_dict.items():
            d[amp_id][ms_id] = dict()
            ms = Microsatellite.objects.get(pk=ms_id)
            for called_alleles, repeat_numbers_dict in genotypes_dict.items():
                dbhist = Histogram.objects.create(
                    sample_reads=dummy_samplereads,
                    microsatellites_version=0,
                    amplicon_id=amp_id,
                    num_reads=sum(repeat_numbers_dict.values())
                )
                # So our objects don't have "special" objects in fields
                dbhist = Histogram.objects.get(pk=dbhist.pk)
                for repeat_number, reads_number in repeat_numbers_dict.items():
                    f1 = get_unique_path("fastq")
                    touch(f1)
                    f2 = get_unique_path("fastq")
                    touch(f2)
                    fm = get_unique_path("fastq")
                    touch(fm)
                    mhg, created = MicrosatelliteHistogramGenotype.objects.get_or_create(
                        microsatellite=ms,
                        repeat_number=repeat_number,
                    )
                    her = HistogramEntryReads.objects.create(
                        histogram=dbhist,
                        microsatellite_genotypes=MicrosatelliteHistogramGenotypeSet.get_for_msgs([mhg]),
                        snp_genotypes=snp_histogram_genotypes,
                        num_reads=reads_number,
                        fastq1=f1,
                        fastq2=f2,
                        fastqm=fm,
                    )
                d[amp_id][ms_id][called_alleles] = dbhist
    yield d
    Histogram.objects.all().delete()
