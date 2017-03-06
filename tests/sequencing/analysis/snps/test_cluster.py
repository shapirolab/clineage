import pytest
import os
from Bio import SeqIO
from frogress import bar
from sequencing.analysis.snps.snp_calling import run_parallel, snp_table

from distributed.utils_test import cluster
from distributed import as_completed, Executor
from distributed.client import Future
from sequencing.analysis.snps.models import ReadsAlignment, VCFReads, SNPReads
from sequencing.analysis.models import AmpliconCollectionBWAIndex
from sequencing.analysis.models_common import BWAIndexMixin
from targeted_enrichment.amplicons.models import AmpliconCollection

import django


@pytest.yield_fixture(scope="session")
def executor():
    with cluster(4) as (d, workers):
        e = Executor(("127.0.0.1", d["port"]))
        yield e
        e.shutdown()


@pytest.mark.django_db
def test_run_parallel(executor,
                      demultiplexing,
                      snp_reads_d,
                      requires_snps,
                      requires_snp_targets,
                      requires_none_genotypes,
                      temp_storage):

    dm = demultiplexing.samplereads_set.all()[0]
    amplicon_collection = dm.library.subclass.panel.amplicon_collection
    parallel_list = run_parallel(executor, demultiplexing.samplereads_set.all(),
                                 amplicon_collection, min_cover=1)

    read_alignment = next(parallel_list)
    vcf_list = next(parallel_list)
    snp_list = next(parallel_list)
    snp_files_list = []
    for snp_future in snp_list:
        snp_f = snp_future.result()
        snp_files_list.append(snp_f)

    snp_dict, rows = snp_table(snp_files_list, min_cover=1)

    for snp_future in snp_list:
        snp_f = snp_future.result()
        snp_obj = snp_f.extract_snp_file()
        cell_id = snp_f.vcf_read.reads_alignment.sample_read.id
        assert snp_dict[('6', '112')][str(cell_id)]['base'] == 'G'

        for chrom, loc in snp_obj:
            assert snp_dict[(chrom, loc)][str(cell_id)]['base'] == snp_obj[(chrom, loc)]['base']

    for Model in [
        SNPReads,
        VCFReads,
        ReadsAlignment,
        AmpliconCollectionBWAIndex,
    ]:
        Model.objects.all().delete()
