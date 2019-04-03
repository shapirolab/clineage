from sequencing.analysis.models import SampleReads, AdamHistogram, FullMSVHistogram
from sequencing.runs.models import NGSRun, Demultiplexing
from sampling.models import Individual, SampleComposition
from lib_prep.workflows.models import MagicalOM6BarcodedContent, MagicalPCR1BarcodedContent
from sequencing.calling.queries.mutation_maps import get_cas_dict, transpose_dict
from sequencing.calling.simcor.calling import split_genotypes
from sequencing.calling.queries.order_this import invert_calling_assignments, binned_genotypes_bi

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

def get_multi_run_srs(ind, run_names,single_cell_only=True, group_label_only=True):
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
        skiped = [sr for sr in srs if not sr.cell.custom_group_labeling]
        if skiped:
            print('WARNING: removing {} cells with no group label'.format(len(skiped)))
        if group_label_only:
            srs = [sr for sr in srs if sr.cell.custom_group_labeling]
        # if demux_scheme_name == 'demux_scheme':
        #     srs = [sr for sr in srs if sr.cell.custom_group_labeling in
        #            ['Kera','Bulk']
        #           ]
        multi_run_srs += srs
    return multi_run_srs, histogram_class


def get_classifying_microsatelites(ind, multi_run_srs, calling_scheme, reads_threshold, confidence_threshold,\
                                   histogram_class, max_distance_from_peak, minimal_distance_between_peaks, \
                                   case, filter_ones, filter_single_peak):
    classifying_ms = dict()
    bin_map_by_ms = dict()
    if calling_scheme.allele_number == 1:
        mono_scheme = True
    else:
        mono_scheme = False

    cas_d = get_cas_dict(srs=multi_run_srs, calling_scheme=calling_scheme, confidence_threshold=confidence_threshold,\
                         reads_threshold=reads_threshold, histogram_class=histogram_class)
    cas_d_by_ms = transpose_dict(cas_d)
    for ms in cas_d_by_ms:
        if not mono_scheme and (ms.slice.chromosome.name not in ['X','Y'] or ind.sex == 'F'):
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
        elif mono_scheme and ms.slice.chromosome.name in ['X', 'Y'] and ind.sex == 'M':
            genotypes = [a for ca in cas_d_by_ms[ms].values() for a in ca.genotypes.alleles]
            genotypes_iterable = [genotypes]

        # other possible cases that are illegal
        elif mono_scheme and ms.slice.chromosome.name not in ['X', 'Y']:
            continue
        elif not mono_scheme and ms.slice.chromosome.name in ['X', 'Y'] and ind.sex == 'M':
            continue
        elif mono_scheme and ind.sex == 'F':
            continue
        else:
            raise invalid_calling_situation("invalid_calling_situation. allele number is {}, individual sex is {} and chromosome is {}".format(\
                calling_scheme.allele_number,ind.sex,ms.slice.chromosome.name))

        for genotypes in genotypes_iterable:
            top_classifying_alleles = tuple(sorted(Counter(genotypes).values(), reverse=True)[:2])
            if len(top_classifying_alleles) <= 1:
                continue
            classifying_ms.setdefault(top_classifying_alleles, list()).append(ms)

    return classifying_ms, bin_map_by_ms

