from targeted_enrichment.planning.models import Microsatellite
from frogress import bar
from sequencing.calling.simcor.calling import split_genotypes
from sequencing.calling.simcor.calling import ms_genotypes_population_query_with_amplicon_all
from sequencing.analysis.models import Histogram
import numpy as np


def transpose_dict(d):
    td = dict()
    for k1 in d:
        for k2 in d[k1]:
            td.setdefault(k2, dict())[k1] = d[k1][k2]
    return td


def filter_mutation_map(mutations_dict, pl1=50, pl2=50):
    l1_counts = {kl1: len(mutations_dict[kl1]) for kl1 in mutations_dict}
    trans_mutations_dict = transpose_dict(mutations_dict)
    l2_counts = {kl2: len(trans_mutations_dict[kl2]) for kl2 in trans_mutations_dict}
    l1_lim = np.percentile(list(l1_counts.values()), pl1)
    l2_lim = np.percentile(list(l2_counts.values()), pl2)
    filtered_by_l1 = {kl1: d for kl1, d in mutations_dict.items() if l1_counts[kl1] > l1_lim}
    filtered_by_l1_l2 = transpose_dict({kl2: d for kl2, d in transpose_dict(filtered_by_l1).items() if l2_counts[kl2] > l2_lim})
    return filtered_by_l1_l2


def get_mono_mutations_dict(srs, calling_scheme, confidence_threshold=0.01, reads_threshold=30,
                            ms_repaet_unit='AC', histogram_class=Histogram):
    amp_by_ms = map_amplicons_to_ms(srs)
    d = dict()
    for ms, amp in bar(amp_by_ms.items()):
        if ms.repeat_unit_type != ms_repaet_unit:
            continue
        for ca in ms_genotypes_population_query_with_amplicon_all(ms, amp, srs, calling_scheme,
                                                        confidence=confidence_threshold,
                                                        reads_threshold=reads_threshold,
                                                        histogram_class=histogram_class):
            if ca.confidence > confidence_threshold:
                continue
            d.setdefault(ca.histogram.sample_reads, dict())[ca.microsatellite] = ca.genotypes.allele1
    return d


def map_amplicons_to_ms(srs):
    # Map the amplicons to microsatellites so that we can filter histograms by ms
    panels = set(sr.library.subclass.panel for sr in srs)
    assert len(panels) == 1
    panel = panels.pop()
    amp_by_ms = dict()
    for amp in panel.amplicon_collection.amplicons.all():
        for ms in Microsatellite.objects.filter(slice__contained=amp.slice):
            amp_by_ms[ms] = amp
    return amp_by_ms


def get_bi_mutations_dict(srs, calling_scheme, confidence_threshold=0.01, reads_threshold=30, max_distance_from_peak=3,
                          ms_repaet_unit='AC', histogram_class=Histogram):
    amp_by_ms = map_amplicons_to_ms(srs)
    ms_split_calling_results = dict()
    for ms, amp in bar(amp_by_ms.items()):
        if ms.repeat_unit_type != ms_repaet_unit:
            continue
        ms_split_calling_results[ms] = split_genotypes(ms, srs, amp, calling_scheme,
                                                       max_distance_from_peak=max_distance_from_peak,
                                                       confidence=confidence_threshold,
                                                       reads_threshold=reads_threshold,
                                                       histogram_class=histogram_class)
    return ms_split_calling_results


def invert_biallelic_results(ms_split_calling_results):
    """
    Inverting a nested dictionary from:
        d[sr][ms][bin][allele]
     to:
        d[ms][sr][allele][bin]
    Args:
        ms_split_calling_results: 

    Returns:

    """
    tran_h1 = dict()
    for ms in ms_split_calling_results:
        if ms_split_calling_results[ms] is None:
            tran_h1[ms.id] = None
            continue
        tran_h1[ms] = dict()
        for ca in ms_split_calling_results[ms]:
            inverted_dict = dict()
            for allele, slot in ms_split_calling_results[ms][ca].items():
                inverted_dict.setdefault(slot, []).append(allele)
            tran_h1[ms][ca.histogram.sample_reads] = {k: v[0] for k, v in inverted_dict.items() if len(v) == 1}
    return tran_h1


def add_mono_calling_to_hemizygous_loci(biallelic_dict, mono_dict, hemizygous_chromosomes=frozenset({'X', 'Y'})):
    tran_ms_mono_and_bi = invert_biallelic_results(biallelic_dict)
    #add for each ms,sr the mono allelic calling
    for sr in mono_dict:
        for ms in mono_dict[sr]:
            if ms.slice.chromosome.name not in hemizygous_chromosomes:
                continue
            if ms not in tran_ms_mono_and_bi or tran_ms_mono_and_bi[ms] is None:
                tran_ms_mono_and_bi[ms] = dict()
            tran_ms_mono_and_bi.setdefault(ms, dict()).setdefault(sr, dict())['mono'] = mono_dict[sr][ms]
    return tran_ms_mono_and_bi


def flatten_bi_allelic_binning(tran_ms_mono_and_bi):
    """
    Deobjectify SampleReads and Microsatellite objects and leave their ids
    Flatten biallelic bins to ms-like identities MSID_BINI
    Args:
        tran_ms_mono_and_bi: 

    Returns:

    """
    print_ready = dict()
    for ms in tran_ms_mono_and_bi:
        if tran_ms_mono_and_bi[ms] is None:
            continue
        bins = {k for sr in tran_ms_mono_and_bi[ms] for k in tran_ms_mono_and_bi[ms][sr].keys()}
        if 'mono' in bins:  # override biallic calling with mono allelic calling results across all sample reads
            for sr in tran_ms_mono_and_bi[ms]:
                if 'mono' not in tran_ms_mono_and_bi[ms][sr]:
                    continue
                print_ready.setdefault(sr.id, dict())[ms.id] = tran_ms_mono_and_bi[ms][sr]['mono']
        else:
            bins_to_ms_labels = {bin_key: '{}_{}'.format(ms.id, label_i) for label_i, bin_key in enumerate(sorted(list(bins)))}
            for sr in tran_ms_mono_and_bi[ms]:
                for bin_key in tran_ms_mono_and_bi[ms][sr]:
                    print_ready.setdefault(sr.id, dict())[bins_to_ms_labels[bin_key]] = tran_ms_mono_and_bi[ms][sr][bin_key]
    return print_ready
