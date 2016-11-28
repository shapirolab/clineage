
import itertools
import functools
from distributed import as_completed

from sequencing.analysis.full_msv.full_msv import merge,\
    align_reads_to_ms_variations, separate_reads_by_genotypes


def list_iterator(f):
    @functools.wraps(f)
    def inner(*args, **kwargs):
        return list(f(*args, **kwargs))
    return inner


def close_connection_and(f):
    @functools.wraps(f)
    def inner(*args, **kwargs):
        from django.db import connection
        connection.close()
        return f(*args, **kwargs)
    return inner


def run_parallel(executor, sample_reads, included_reads="M", mss_version=0, ref_padding=50):
    # TODO: set resource.getrlimit(resource.RLIMIT_CORE) to something low for all bowtie2 related jobs
    # *currently in dworker.q
    merged_reads = executor.map(merge, sample_reads, pure=False)
    fmsvas = executor.map(
        close_connection_and(align_reads_to_ms_variations), merged_reads,
        itertools.repeat(ref_padding), itertools.repeat(mss_version), pure=False
    )

    fhers_list = executor.map(list_iterator(close_connection_and(separate_reads_by_genotypes)), fmsvas, pure=False)
    yield merged_reads, fmsvas, fhers_list
