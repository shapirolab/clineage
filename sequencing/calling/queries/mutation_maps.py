from sequencing.calling.models import CalledAlleles
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

