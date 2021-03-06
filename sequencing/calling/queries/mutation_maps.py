from targeted_enrichment.planning.models import Microsatellite
from sequencing.calling.simcor.calling import split_genotypes
from sequencing.calling.models import BestCorrelationCalledAlleles
from sequencing.calling.hist import Histogram as dHistogram
from collections import Counter
from sequencing.calling.simcor.hist_analysis import get_far_apart_highest_peaks, better_get_far_apart_highest_peaks
from sequencing.analysis.models import Histogram
import numpy as np


def transpose_dict(d):
    td = dict()
    for k1 in d:
        for k2 in d[k1]:
            td.setdefault(k2, dict())[k1] = d[k1][k2]
    return td


def get_cas_dict(srs, calling_scheme, confidence_threshold=0.01, reads_threshold=30, histogram_class=Histogram):
    d = dict()
    for ca in BestCorrelationCalledAlleles.objects.filter(
        calling_scheme=calling_scheme,
        histogram__in=histogram_class.objects.filter(
            num_reads__gte=reads_threshold,
            sample_reads__in=srs),
        confidence__lte=confidence_threshold).select_related(
            'histogram__sample_reads', 'microsatellite__slice__chromosome', 'genotypes'):
        d.setdefault(ca.histogram.sample_reads, dict())[ca.microsatellite] = ca
    return d


def filter_mutation_map(mutations_dict, pl1=50, pl2=50):
    l1_counts = {kl1: len(mutations_dict[kl1]) for kl1 in mutations_dict}
    trans_mutations_dict = transpose_dict(mutations_dict)
    l2_counts = {kl2: len(trans_mutations_dict[kl2]) for kl2 in trans_mutations_dict}
    l1_lim = np.percentile(list(l1_counts.values()), pl1)
    l2_lim = np.percentile(list(l2_counts.values()), pl2)
    filtered_by_l1 = {kl1: d for kl1, d in mutations_dict.items() if l1_counts[kl1] >= l1_lim}
    filtered_by_l1_l2 = transpose_dict({kl2: d for kl2, d in transpose_dict(filtered_by_l1).items() if l2_counts[kl2] >= l2_lim})
    return filtered_by_l1_l2


def get_mono_mutations_dict(srs, calling_scheme, confidence_threshold=0.01, reads_threshold=30,
                            histogram_class=Histogram):
    cas_d = get_cas_dict(srs, calling_scheme, confidence_threshold=confidence_threshold,
                 reads_threshold=reads_threshold, histogram_class=histogram_class)
    d = dict()
    for sr in cas_d:
        for ms, ca in cas_d[sr].items():
            d.setdefault(ca.histogram.sample_reads, dict())[ca.microsatellite] = ca.genotypes.allele1
    return d


def merge_mono_mutations_dicts(d1, d2):
    """
    Merge two mono allelic mutation dictionaries
    asserts no overlap between the mss
    Args:
        d1: 
        d2: 

    Returns:

    """
    merged = dict()
    for sr in d2.keys() | d1.keys():  # initialize l2 dictionaries
        merged[sr] = dict()
    for sr in d2.keys() - d1.keys():
        merged[sr].update(d2[sr])
    for sr in d1.keys() - d2.keys():
        merged[sr].update(d1[sr])
    for sr in d1.keys() & d2.keys():  # l1 shared keys undergo l2 check for mutual exclusivity
        assert d1[sr].keys() & d2[sr].keys() == set()
        merged[sr].update(d1[sr])
        merged[sr].update(d2[sr])
    return merged


def map_amplicons_to_ms(srs):
    # Map the amplicons to microsatellites so that we can filter histograms by ms
    panels = set(sr.library.subclass.panel for sr in srs)
    assert len(panels) == 1
    panel = panels.pop()
    amps_by_ms = dict()
    for amp in panel.amplicon_collection.amplicons.all():
        for ms in Microsatellite.objects.filter(slice__contained=amp.slice):
            amps_by_ms.setdefault(ms, list()).append(amp)
    return amps_by_ms


def get_bi_mutations_dict(srs, calling_scheme, confidence_threshold=0.01, reads_threshold=30, max_distance_from_peak=3,
                          histogram_class=Histogram, case=1, filter_ones=False):
    cas_d = get_cas_dict(srs, calling_scheme, confidence_threshold=confidence_threshold,
                         reads_threshold=reads_threshold, histogram_class=histogram_class)
    cas_d_by_ms = transpose_dict(cas_d)
    ms_split_calling_results = dict()
    for ms in cas_d_by_ms:
        calling_assignments = split_genotypes(cas_d_by_ms[ms].values(),
                        max_distance_from_peak=max_distance_from_peak, case=case,
                        filter_ones=filter_ones)
        if calling_assignments is None:
            continue
        ms_split_calling_results[ms] = calling_assignments
    return ms_split_calling_results


def invert_biallelic_results(ms_split_calling_results, strict=True):
    """
    Inverting a nested dictionary from:
        d[sr][ms][bin][allele]
     to:
        d[ms][sr][allele][bin]
    Args:
        ms_split_calling_results: 
        strict: if True, assume one allele per bin, drop other cases
                if False, use the average

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
            if strict:
                tran_h1[ms][ca.histogram.sample_reads] = {k: v[0] for k, v in inverted_dict.items() if len(v) == 1}
            else:
                tran_h1[ms][ca.histogram.sample_reads] = {k: sum(v)/len(v) for k, v in inverted_dict.items()}
    return tran_h1


def add_mono_calling_to_hemizygous_loci(biallelic_dict, mono_dict,
                                        hemizygous_chromosomes=frozenset({'X', 'Y'}), strict_alleles=True):
    tran_ms_mono_and_bi = invert_biallelic_results(biallelic_dict, strict=strict_alleles)
    hemizygous_mss = set(Microsatellite.objects.filter(slice__chromosome__name__in=hemizygous_chromosomes))
    #add for each ms,sr the mono allelic calling
    for sr in mono_dict:
        for ms in mono_dict[sr]:
            if ms not in hemizygous_mss:
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


def is_population_mono_allelic(alleles, minimal_distance_between_peaks=2, min_prop=0.05, case=1):
    h = dHistogram(Counter(alleles))
    if case == 1:
        return len(get_far_apart_highest_peaks(
                h,
                allele_number=2,
                minimal_distance_between_peaks=minimal_distance_between_peaks,
                min_prop=min_prop)) == 1
    elif case == 2:
        return len(better_get_far_apart_highest_peaks(
            h,
            minimal_distance_between_peaks=1,  # override with strict value
            min_prop=min_prop)) == 1


def filter_bipartition_loci(full_td, case=1):
    btd = dict()
    tbtd = transpose_dict(full_td)
    for loc in tbtd:
        if is_population_mono_allelic(tbtd[loc].values(), case=case):
            for c in tbtd[loc]:
                btd.setdefault(c, dict())[loc] = full_td[c][loc]
    return btd


from decimal import Decimal
def query_rounded_proportional_mutation_map(ind, multi_run_srs, schema_by_repeat_type, histogram_class):
    mutations_dict = dict()
    for repeat_type, (calling_schema, rt, ct) in schema_by_repeat_type.items():
        cas_d = get_cas_dict(multi_run_srs, calling_schema, confidence_threshold=ct,
                         reads_threshold=rt, histogram_class=histogram_class)
        cas_d_by_ms = transpose_dict(cas_d)
        for ms in cas_d_by_ms:
            if repeat_type in ['ac_b', 'ac_bb'] and (ms.slice.chromosome.name not in ['X', 'Y'] or ind.sex == 'F'):
                for sr in cas_d_by_ms[ms]:
                    alleles = cas_d_by_ms[ms][sr].genotypes.alleles
                    if len(alleles) == 1:
                        prop = Decimal('1.0')
                    elif len(alleles) == 2:
                        prop = Decimal('0.5')
                    else:
                        raise ValueError(len(alleles))
                    mutations_dict.setdefault(ms,dict())[sr]=frozenset(
                        {(int(a), prop) for a in alleles})
            elif ms.slice.chromosome.name in ['X', 'Y'] and repeat_type not in ['ac_b', 'ac_bb'] and ind.sex == 'M':
                for sr in cas_d_by_ms[ms]:
                    alleles = cas_d_by_ms[ms][sr].genotypes.alleles
                    assert len(alleles) == 1
                    mutations_dict.setdefault(ms,dict())[sr]=frozenset(
                        {(int(a), Decimal('1.0')) for a in alleles})
#                 genotypes = [a for ca in cas_d_by_ms[ms].values() for a in ca.genotypes.alleles]
#                 genotypes_iterable = [genotypes]
            elif calling_schema.allele_number == 1 and ms.slice.chromosome.name not in ['X','Y']:
                continue
            elif calling_schema.allele_number == 2 and ms.slice.chromosome.name in ['X','Y'] and ind.sex == 'M':
                continue
            elif calling_schema.allele_number == 1 and ind.sex == 'F':
                continue
            else:
                # sanity check
                print("ind {} did not pass any if condition with scheme {} and ms {}".format(ind.name, repeat_type, ms.id))
                continue
    return mutations_dict


def get_root_genotypes(d_by_ms, multi_run_srs, p_mono=0.7, order=None):
    # p_mono aims to classify between true mono-allelic loci and false ones
    # A good test for that is the number of resulting mono and bi loci from the X chromosome of male VS female
    if order is None:
        root_apx_srs = {sr:sr.barcoded_content.subclass.amplified_content.cell.is_root_approximating for sr in multi_run_srs}
        order = [[sr for sr in multi_run_srs if root_apx_srs[sr]]]
    root_genotypes = dict()
    for ms in d_by_ms:
        root_genotype=None
        for srs in order:
            if ms in root_genotypes:
                continue
            bi_alleles = []
            mono_alleles = []
            for sr in set(srs)&d_by_ms[ms].keys():
                alleles = [a for a,p in d_by_ms[ms][sr]]
                if len(alleles) < 2:
                    mono_alleles.append(sorted(alleles)[0])
                    continue
                bi_alleles.append(tuple(sorted(alleles)))
            if bi_alleles:
                bi_c = Counter(bi_alleles)
                root_genotype = frozenset(
                            {(int(a), Decimal('0.5')) for a in bi_c.most_common()[0][0]})
            if mono_alleles:
                mono_c = Counter(mono_alleles)
                root_genotype = frozenset(
                            {(mono_c.most_common()[0][0], Decimal('1.0'))})
            if bi_alleles and mono_alleles:
                if mono_c.most_common()[0][1]/(bi_c.most_common()[0][1]+mono_c.most_common()[0][1]) > p_mono:
                    root_genotype = frozenset(
                            {(mono_c.most_common()[0][0], Decimal('1.0'))})
                else:
                    root_genotype = frozenset(
                            {(int(a), Decimal('0.5')) for a in bi_c.most_common()[0][0]})
            if root_genotype is None:
                continue
            root_genotypes[ms] = root_genotype
    return root_genotypes