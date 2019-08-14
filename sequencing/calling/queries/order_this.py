from sampling.models import SampleComposition
from sequencing.calling.hist import Histogram as dHistogram
from sequencing.calling.simcor.hist_analysis import better_get_far_apart_highest_peaks
from collections import Counter
from sampling.models import Individual
from lib_prep.workflows.models import MagicalOM6BarcodedContent
from sequencing.calling.simcor.calling import split_genotypes
from sequencing.calling.models import BestCorrelationCalledAlleles
from sequencing.runs.models import NGSRun
from sequencing.analysis.models import SampleReads
from sequencing.analysis.full_msv.models import FullMSVHistogram
from sequencing.calling.queries.mutation_maps import get_cas_dict, transpose_dict


def query_individual_and_sample_reads(ind_name, run_name, demux_class):
    sc = SampleComposition.objects.get(name='Single Cell')
    ind = Individual.objects.get(name=ind_name)
    run = NGSRun.objects.get(name=run_name)
    try:
        demux = demux_class.objects.get(ngs_run=run)
    except demux_class.DoesNotExist:
        demux = demux_class.objects.filter(ngs_run=run).order_by('-pk')[0]  # latest demux
    bcs = MagicalOM6BarcodedContent.objects.filter(content__cell__individual=ind)
    srs = SampleReads.objects.filter(demux=demux).filter(barcoded_content__in=bcs)
    srs = [sr for sr in srs if sr.barcoded_content.subclass.amplified_content.cell.composition == sc]
    return ind, srs


def split_heterozygous_identities(ms_cas_d):
    calling_assignments = split_genotypes(ms_cas_d.values(),
                    max_distance_from_peak=2,
                    case=2,
                    filter_ones=True,
                    filter_single_peak=True,)
    if calling_assignments is None:
        return None, None
    ica = invert_calling_assignments(calling_assignments)
    genotypes_iterable = binned_genotypes_bi(ica)
    return ica, genotypes_iterable


def get_peaks(genotypes):
    h = dHistogram(
                {k: v for k, v in Counter(genotypes).items() if v > 1}
            )
    if len(h.keys()) == 0:
        return None
    peaks = better_get_far_apart_highest_peaks(
        h,
        minimal_distance_between_peaks=1,
        min_prop=0.05)
    return peaks


def invert_calling_assignments(calling_assignments, average=False):
    ica = dict()
    for ca in calling_assignments:
        inverted_dict = dict()
        for allele, slot in calling_assignments[ca].items():
            inverted_dict.setdefault(slot, []).append(allele)
        if average:
            ica[ca.histogram.sample_reads] = {k: sum(v)//len(v) for k, v in inverted_dict.items()}
        else:
            ica[ca.histogram.sample_reads] = {k: v[0] for k, v in inverted_dict.items() if len(v) == 1}
    return ica


def binned_genotypes_bi(ica):
    bins = {k for sr in ica for k in ica[sr].keys()}
    for bin_key in bins:
        genotypes = []
        for sr in ica:
            if bin_key in ica[sr]:
                genotypes.append(ica[sr][bin_key])
        yield genotypes


def get_viable_mss(srs, repeat_type, schema_by_repeat_type, drop_multiple_peaks=True, drop_unclassifying_alleles=False, debug_prints=False):
    calling_schema, rt, ct = schema_by_repeat_type[repeat_type]
    cas_d = get_cas_dict(srs, calling_schema, confidence_threshold=ct,
                     reads_threshold=rt, histogram_class=FullMSVHistogram)  # bulk query called alleles
    cas_d_by_ms = transpose_dict(cas_d)
    bin_map_by_ms = dict()
    classifying_ms_by_major_pair = dict()
    for ms in cas_d_by_ms:
        if ms.slice.chromosome.name not in ['X', 'Y'] and repeat_type in ['ac_b']:
            ica, genotypes_iterable = split_heterozygous_identities(cas_d_by_ms[ms])
            bin_map_by_ms[ms] = ica
        elif ms.slice.chromosome.name in ['X', 'Y'] and repeat_type not in ['ac_b']:
            genotypes = [a for ca in cas_d_by_ms[ms].values() for a in ca.genotypes.alleles]
            genotypes_iterable = [genotypes]
        else:
            continue
        if genotypes_iterable is None:
            continue  # haplotyping failure
        for genotypes in genotypes_iterable:
            peaks = get_peaks(genotypes)
            if peaks is not None and len(peaks) != 1 and drop_multiple_peaks:
                if max(peaks) - min(peaks) > 3:
                    if debug_prints:
                        print(
                            ms,
                            'multiple distanct peaks',
                            repeat_type,
                            Counter(genotypes),
                            Counter([a for ca in cas_d_by_ms[ms].values() for a in ca.genotypes.alleles]))
                    continue  # Drop loci with multiple distinct peaks
            top_classifying_alleles = tuple(sorted(Counter(genotypes).values())[:2])
            if len(top_classifying_alleles) <= 1 and drop_unclassifying_alleles:
                continue  # Drop loci that show the same signal across all cells
            classifying_ms_by_major_pair.setdefault(top_classifying_alleles, list()).append(ms)
    return classifying_ms_by_major_pair, bin_map_by_ms


def get_mutation_map(srs, classifying_ms_by_type, bin_map_by_ms, schema_by_repeat_type):
    d = dict()
    for repeat_type in schema_by_repeat_type:
        if repeat_type == 'ac_b':
            for ms in bin_map_by_ms.keys() & set(classifying_ms_by_type[repeat_type]):
                for sr in bin_map_by_ms[ms].keys() & set(srs):
                    for bin_key, allele in bin_map_by_ms[ms][sr].items():
                        d.setdefault(sr, dict())[(ms, bin_key)] = allele
        else:
            (schema, rt, ct) = schema_by_repeat_type[repeat_type]
            for ca in BestCorrelationCalledAlleles.objects.filter(
                calling_scheme=schema,
                microsatellite__in=classifying_ms_by_type[repeat_type],
                histogram__in=FullMSVHistogram.objects.filter(
                    num_reads__gte=rt,
                    sample_reads__in=srs),
                confidence__lte=ct).select_related(
                    'histogram__sample_reads', 'microsatellite__slice__chromosome', 'genotypes'):
                d.setdefault(ca.histogram.sample_reads, dict())[ca.microsatellite] = ca.genotypes.allele1
    return d
