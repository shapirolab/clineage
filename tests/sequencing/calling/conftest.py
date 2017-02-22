import pytest
import os
import itertools
from Bio import SeqIO

from sequencing.analysis.full_msv.models import SampleReads, FullMSVMergedReads, \
    FullMSVAssignment, FullMSVHistogram, FullMSVariations, Histogram
from sequencing.analysis.models_common import PearOutputMixin, SNPHistogramGenotypeSet, \
    MicrosatelliteHistogramGenotypeSet, SNPHistogramGenotype, MicrosatelliteHistogramGenotype, \
    HistogramEntryReads
from misc.utils import get_unique_path
from sequencing.analysis.full_msv.full_msv import get_full_ms_variations

from tests.sequencing.runs.conftest import *
from tests.sequencing.analysis.full_msv.full_msv_amp1_3_4 import VARS_134
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.amplicons.conftest import *
from sequencing.calling.models import CallingScheme
from tests.sequencing.calling.reads_dict import MS_HISTOGRAMS_DICT
from tests.flat_dict import FlatDict


def touch(fp):
    with open(fp, 'wb') as f:
        pass


@pytest.fixture(scope="session")
def histograms_fd():
    return FlatDict(MS_HISTOGRAMS_DICT)


@pytest.fixture()
def requires_none_genotypes(request, transactional_db):
    MicrosatelliteHistogramGenotype.objects.get_or_create(microsatellite=None,
        defaults=dict(repeat_number=1))
    SNPHistogramGenotype.objects.get_or_create(snp=None, defaults=dict(base=""))


@pytest.fixture()
def dummy_samplereads(demultiplexing):
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
    yield sr
    sr.delete()


def histograms_and_calling_solutions(dummy_samplereads, histograms_fd, requires_microsatellites, requires_none_genotypes):
    none_snp_genotype = SNPHistogramGenotype.objects.get(snp=None)
    snp_histogram_genotypes, c = SNPHistogramGenotypeSet.objects.get_or_create(
        **{fn: none_snp_genotype for fn in SNPHistogramGenotypeSet.genotype_field_names()})
    for amp_id, ms_dict in histograms_fd.items():
        for ms_id, genotypes_dict in ms_dict.items():
            ms = Microsatellite.objects.get(pk=ms_id)
            for called_alleles, repeat_numbers_dict in genotypes_dict.items():
                dbhist = Histogram.objects.create(
                    sample_reads=dummy_samplereads,
                    microsatellites_version=0,
                    amplicon=amp_id,
                    num_reads=sum(ms_dict.values())
                )
                for repeat_number, reads_number in repeat_numbers_dict.items():
                    f1 = get_unique_path("fastq")
                    touch(f1)
                    f2 = get_unique_path("fastq")
                    touch(f2)
                    fm = get_unique_path("fastq")
                    touch(fm)
                    her = HistogramEntryReads.objects.create(
                        histogram=dbhist,
                        microsatellite_genotypes=MicrosatelliteHistogramGenotypeSet.get_for_msgs(
                            [MicrosatelliteHistogramGenotype.objects.get_or_create(
                                microsatellite=ms,
                                repeat_number=repeat_number,
                            )]
                        ),
                        snp_genotypes=snp_histogram_genotypes,
                        num_reads=reads_number,
                        fastq1=f1,
                        fastq2=f2,
                        fastqm=fm,
                    )
                yield called_alleles, dbhist, ms
                dbhist.delete()
