from targeted_enrichment.planning.models import Microsatellite
from frogress import bar
from sequencing.calling.simcor.calling import split_genotypes
from sequencing.calling.models import CalledAlleles
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


def get_mono_mutations_dict(srs, calling_scheme, confidence_threshold=0.01, reads_threshold=30):
    d = dict()
    for sr in srs:
        d[sr] = dict()
        for ca in CalledAlleles.objects.filter(
                calling_scheme=calling_scheme,
                histogram__sample_reads=sr,
                histogram__num_reads__gte=reads_threshold).select_subclasses():
            if ca.confidence > confidence_threshold:
                continue
            d[sr][ca.microsatellite] = ca.genotypes.allele1
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
