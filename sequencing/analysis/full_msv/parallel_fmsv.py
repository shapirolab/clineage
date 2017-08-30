
import itertools
import functools
from distributed import as_completed

from sequencing.analysis.full_msv.full_msv import merge,\
    align_reads_to_ms_variations, separate_reads_by_genotypes, \
    split_merged_reads_as_list, align_reads_to_ms_variations_part,\
    merge_fmsva_parts


def double_map(executor, func, future_lists, *params):
    for f in as_completed(future_lists):
        try:
            l = f.result()
        except:
            yield []
        else:
            yield executor.map(func, l, *[itertools.repeat(p) for p in params], pure=False)


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


def run_parallel(executor, sample_reads, included_reads="M", mss_version=0, ref_padding=50, amp_col_size=15000):
    # TODO: set resource.getrlimit(resource.RLIMIT_CORE) to something low for all bowtie2 related jobs
    # *currently in dworker.q
    merged_reads = executor.map(merge, sample_reads, pure=False)
    fmsvas = executor.map(
        close_connection_and(align_reads_to_ms_variations), merged_reads,
        itertools.repeat(ref_padding), itertools.repeat(mss_version), itertools.repeat(amp_col_size), pure=False
    )

    fhers_list = executor.map(list_iterator(close_connection_and(separate_reads_by_genotypes)), fmsvas, pure=False)
    yield merged_reads, fmsvas, fhers_list


def run_parallel_split_alignments(executor, sample_reads, included_reads="M", reads_chunk_size=10**5, mss_version=0, ref_padding=50):
    # TODO: set resource.getrlimit(resource.RLIMIT_CORE) to something low for all bowtie2 related jobs
    # *currently in dworker.q
    merged_reads = executor.map(merge, sample_reads, pure=False)
    yield merged_reads
    fmsv_merged_reads_parts_lists = executor.map(
        close_connection_and(split_merged_reads_as_list),
        merged_reads,
        itertools.repeat(reads_chunk_size),
        itertools.repeat(included_reads),
        pure=False
    )
    yield fmsv_merged_reads_parts_lists
    fmsva_parts_lists = double_map(executor, close_connection_and(align_reads_to_ms_variations_part),
                                 fmsv_merged_reads_parts_lists, ref_padding, mss_version)
    yield fmsva_parts_lists
    merged_fmsvas = executor.map(
        close_connection_and(merge_fmsva_parts),
        fmsva_parts_lists,
        itertools.repeat(reads_chunk_size),
        itertools.repeat(included_reads)
    )
    yield merged_fmsvas
    fhers_list = executor.map(list_iterator(close_connection_and(separate_reads_by_genotypes)), merged_fmsvas, pure=False)
    yield fhers_list
