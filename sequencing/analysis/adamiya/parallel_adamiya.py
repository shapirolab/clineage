
import itertools
import functools
from distributed import as_completed

from sequencing.analysis.adamiya.adamiya import merge, create_reads_index, \
    align_primers_to_reads, separate_reads_by_amplicons,\
    align_reads_to_ms_variations, separate_reads_by_genotypes


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


def run_parallel(executor, sample_reads, included_reads="F", mss_version=0, read_padding=5, ref_padding=50):
    # TODO: set resource.getrlimit(resource.RLIMIT_CORE) to something low for all bowtie2 related jobs
    # *currently in dworker.q
    merged_reads = executor.map(merge, sample_reads, pure=False)
    reads_indices = executor.map(create_reads_index, merged_reads,
        itertools.repeat(included_reads), itertools.repeat(read_padding), pure=False)
    adam_margin_assignments = executor.map(close_connection_and(align_primers_to_reads),
        reads_indices, pure=False)
    adam_amplicon_reads_lists = executor.map(
        list_iterator(close_connection_and(separate_reads_by_amplicons)),
        adam_margin_assignments, pure=False)
    yield merged_reads, reads_indices, adam_margin_assignments, adam_amplicon_reads_lists
    adam_histograms = double_map(executor, close_connection_and(align_reads_to_ms_variations),
        adam_amplicon_reads_lists, ref_padding, mss_version)
    for fs in adam_histograms:
        fs2 = executor.map(list_iterator(close_connection_and(separate_reads_by_genotypes)), fs, pure=False)
        yield fs, fs2
