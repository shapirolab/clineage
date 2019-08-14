from sequencing.analysis.models import SampleReads, AdamHistogram, FullMSVHistogram
from sequencing.runs.models import NGSRun, Demultiplexing
from sampling.models import Individual, SampleComposition
from lib_prep.workflows.models import MagicalOM6BarcodedContent, MagicalPCR1BarcodedContent
from sequencing.calling.queries.mutation_maps import get_cas_dict, transpose_dict
from sequencing.calling.simcor.calling import split_genotypes
from sequencing.calling.queries.order_this import invert_calling_assignments, binned_genotypes_bi
from sequencing.calling.models import BestCorrelationCalledAlleles
from collections import Counter


class invalid_calling_situation(Exception):
    pass


def get_bcs_and_histogram_type(ind, run_name, aa_runs=('nsr4', 'nsr5')):
    if run_name in aa_runs:
        bcs = MagicalPCR1BarcodedContent.objects.filter(content__cell__individual=ind)
        histogram_class = AdamHistogram
    else:
        bcs = MagicalOM6BarcodedContent.objects.filter(content__cell__individual=ind)
        histogram_class = FullMSVHistogram
    return bcs, histogram_class


aa_runs=frozenset({'nsr4', 'nsr5'})  # generally constant unless future aa runs are added


def get_multi_run_srs(ind, run_names, single_cell_only=False, group_label_only=True):
    sc = SampleComposition.objects.get(name='Single Cell')
    multi_run_srs = []
    runs_only = {run_name for run_name, demux_scheme_name in run_names}
    assert runs_only & aa_runs != {} or runs_only - aa_runs == runs_only # same individual can't have runs analyzed differently
    for run_name, demux_scheme_name in run_names:
        run = NGSRun.objects.get(name=run_name)
        demux = Demultiplexing.objects.get(ngs_run=run, demux_scheme__name=demux_scheme_name)
        bcs, histogram_class = get_bcs_and_histogram_type(ind, run_name, aa_runs=aa_runs)
        srs = SampleReads.objects.filter(demux=demux).filter(barcoded_content__in=bcs)
        if single_cell_only:
            srs = [sr for sr in srs if sr.cell.composition == sc]
        skipped = [sr for sr in srs if not sr.cell.custom_group_labeling]
        if group_label_only and len(skipped):
            print('WARNING: removing {} cells with no group label'.format(len(skipped)))
            srs = [sr for sr in srs if sr.cell.custom_group_labeling]
        multi_run_srs += srs
    return multi_run_srs, histogram_class


def cas_d_by_ms_by_sex(
        sex, multi_run_srs, calling_scheme, reads_threshold=10, confidence_threshold=0.02,
        histogram_class=FullMSVHistogram):

    if calling_scheme.allele_number == 1:
        is_mono_scheme = True
    else:
        is_mono_scheme = False

    if sex == 'F' and is_mono_scheme:
        return dict(), is_mono_scheme  # empty dictionaries since female is bi-allelic

    cas_d = get_cas_dict(srs=multi_run_srs, calling_scheme=calling_scheme, confidence_threshold=confidence_threshold,\
                         reads_threshold=reads_threshold, histogram_class=histogram_class)
    cas_d_by_ms = transpose_dict(cas_d)
    return cas_d_by_ms, is_mono_scheme


def filter_classifying_microsatelites(cas_d_by_ms, is_mono_scheme, sex, max_distance_from_peak=2,
                                      minimal_distance_between_peaks=3, case=2, filter_ones=False,
                                      filter_single_peak=False):

    classifying_ms = dict()
    bin_map_by_ms = dict()
    for ms in cas_d_by_ms:
        if not is_mono_scheme and (ms.slice.chromosome.name not in ['X', 'Y'] or sex == 'F'):
            calling_assignments = split_genotypes(cas_d_by_ms[ms].values(),
                                                  max_distance_from_peak=max_distance_from_peak,
                                                  minimal_distance_between_peaks=minimal_distance_between_peaks,
                                                  case=case,
                                                  filter_ones=filter_ones,
                                                  filter_single_peak=filter_single_peak, )

            if calling_assignments is None:
                continue
            ica = invert_calling_assignments(calling_assignments)
            bin_map_by_ms[ms] = ica
            genotypes_iterable = binned_genotypes_bi(ica)
        # mono case
        elif is_mono_scheme and ms.slice.chromosome.name in ['X', 'Y'] and sex == 'M':
            genotypes = [a for ca in cas_d_by_ms[ms].values() for a in ca.genotypes.alleles]
            genotypes_iterable = [genotypes]

        # other possible cases that are illegal
        elif is_mono_scheme and ms.slice.chromosome.name not in ['X', 'Y']:
            continue
        elif not is_mono_scheme and ms.slice.chromosome.name in ['X', 'Y'] and sex == 'M':
            continue
        elif is_mono_scheme and sex == 'F':
            continue
        else:
            raise invalid_calling_situation("invalid_calling_situation. individual sex is {} and chromosome is {}".format(\
                sex, ms.slice.chromosome.name))

        for genotypes in genotypes_iterable:
            top_classifying_alleles = tuple(sorted(Counter(genotypes).values(), reverse=True)[:2])
            if len(top_classifying_alleles) <= 1:
                continue
            classifying_ms.setdefault(top_classifying_alleles, list()).append(ms)

    return classifying_ms, bin_map_by_ms


def cell_labeler(sr):
    return 'ID{}'.format(sr.id)


def query_calling_dict_bi(ind, run_names, repeat_type, bin_map_by_ind_by_ms, classifying_ms_by_ind_by_type, multi_run_srs):
    d = dict()
    for ms in bin_map_by_ind_by_ms[(ind, run_names)].keys() & set(
            classifying_ms_by_ind_by_type[(ind, run_names)][repeat_type]):
        for sr in bin_map_by_ind_by_ms[(ind, run_names)][ms].keys() & set(multi_run_srs):
            sr_label = cell_labeler(sr)
            for bin_key, allele in bin_map_by_ind_by_ms[(ind, run_names)][ms][sr].items():
                d.setdefault(sr_label, dict())['{}_{}_{}'.format(
                    ms.repeat_unit_type,
                    ms.id,
                    bin_key
                )] = allele
    return d


def query_calling_dict_mono(ind, run_names, repeat_type, classifying_ms_by_ind_by_type, multi_run_srs,
                            schema_by_repeat_type, histogram_class=FullMSVHistogram):
    d = dict()
    if ind.sex != 'F' or repeat_type != 'ac':  # don't overwrite biX with monoX
        (schema, rt, ct) = schema_by_repeat_type[repeat_type]
    else:
        return 
    # query cas with enough reads and confidence
    # populate dictionary with sr label as main key, ms label as sub key and allele as value
    for ca in BestCorrelationCalledAlleles.objects.filter(
            calling_scheme=schema,
            microsatellite__in=classifying_ms_by_ind_by_type[(ind, run_names)][repeat_type],
            histogram__in=histogram_class.objects.filter(
                num_reads__gte=rt,
                sample_reads__in=multi_run_srs),
            confidence__lte=ct).select_related(
                'histogram__sample_reads', 'microsatellite__slice__chromosome', 'genotypes'):
        sr_label = cell_labeler(ca.histogram.sample_reads)
        d.setdefault(sr_label, dict())['{}_{}'.format(
            ca.microsatellite.repeat_unit_type,
            ca.microsatellite.id
        )] = ca.genotypes.allele1
    return d


def filter_by_minimal_group_size(classifying_ms_by_ind_by_type_by_grouping, group_size=2):
    classifying_ms_by_ind_by_type = dict()
    for ind, run_names in classifying_ms_by_ind_by_type_by_grouping:
        classifying_ms_by_ind_by_type[(ind, run_names)] = dict()
        for repeat_type in classifying_ms_by_ind_by_type_by_grouping[(ind, run_names)]:
            classifying_ms_by_ind_by_type[(ind, run_names)][repeat_type] = list()
            for tca in classifying_ms_by_ind_by_type_by_grouping[(ind, run_names)][repeat_type]:
                if group_size == 0:
                    classifying_ms_by_ind_by_type[(ind, run_names)][repeat_type] += \
                    classifying_ms_by_ind_by_type_by_grouping[(ind, run_names)][repeat_type][tca]
                    continue
                if len(tca) > 1 and min(tca) >= group_size:
                    classifying_ms_by_ind_by_type[(ind, run_names)][repeat_type] += \
                    classifying_ms_by_ind_by_type_by_grouping[(ind, run_names)][repeat_type][tca]
    return classifying_ms_by_ind_by_type
