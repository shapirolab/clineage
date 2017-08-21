from sequencing.calling.models import HighestPeaksMonoSimCorSchemeModel, \
    HighestPeaksProximityRatioFilteredBiSimCorSchemeModel
from sequencing.calling.queries.mutation_maps import get_mono_mutations_dict, get_bi_mutations_dict
from sequencing.analysis.models import FullMSVHistogram
from sequencing.calling.queries.mutation_maps import merge_mono_mutations_dicts, add_mono_calling_to_hemizygous_loci, \
    flatten_bi_allelic_binning, filter_bipartition_loci, transpose_dict
from targeted_enrichment.planning.models import Microsatellite
from sequencing.calling.queries.formatters import textify_keys_in_mutations_dict
from sequencing.phylo.matlab_wrappers import sr_label_func, ms_label_func


def get_mutation_maps_a_mono(srs, reads_threshold=30, confidence=0.01,
                             schema_class=HighestPeaksMonoSimCorSchemeModel, histogram_class=FullMSVHistogram):
    a_schema_mono = schema_class.objects.get(name='PolyA markov mono',)
    mono_a_dict = get_mono_mutations_dict(
        srs,
        a_schema_mono,
        confidence_threshold=confidence,
        reads_threshold=reads_threshold,
        histogram_class=histogram_class)
    return mono_a_dict


def get_mutation_maps_g_mono(srs, reads_threshold=30, confidence=0.01,
                             schema_class=HighestPeaksMonoSimCorSchemeModel, histogram_class=FullMSVHistogram):
    g_schema_mono = schema_class.objects.get(name='PolyG markov mono',)
    mono_g_dict = get_mono_mutations_dict(
        srs,
        g_schema_mono,
        confidence_threshold=confidence,
        reads_threshold=reads_threshold,
        histogram_class=histogram_class)
    return mono_g_dict


def get_mutation_maps_ac_mono(srs, reads_threshold=30, confidence=0.01,
                              schema_class=HighestPeaksMonoSimCorSchemeModel, histogram_class=FullMSVHistogram):
    ac_schema_mono_extended = schema_class.objects.get(name='AC markov mono extended',)
    mono_ac_dict = get_mono_mutations_dict(
        srs,
        ac_schema_mono_extended,
        confidence_threshold=confidence,
        reads_threshold=reads_threshold,
        histogram_class=histogram_class)
    return mono_ac_dict


def get_mutation_maps_ac_bi(srs, reads_threshold=30, confidence=0.01,
                            schema_class=HighestPeaksProximityRatioFilteredBiSimCorSchemeModel,
                            histogram_class=FullMSVHistogram):
    ac_schema_bi_strict_extended_fast = schema_class.objects.get(name='AC markov HP-PRF strinct extended fast',)
    biallelic_ac_dict = get_bi_mutations_dict(
        srs,
        ac_schema_bi_strict_extended_fast,
        confidence_threshold=confidence,
        reads_threshold=reads_threshold,
        max_distance_from_peak=3,
        histogram_class=histogram_class)
    return biallelic_ac_dict


def combine_cases(mono_a, mono_g, mono_ac, bi_ac):
    if mono_a is not None:
        mono_ac_a = merge_mono_mutations_dicts(mono_ac, mono_a)
    else:
        mono_ac_a = mono_ac
    if mono_g is not None:
        full_mono = merge_mono_mutations_dicts(mono_ac_a, mono_g)
    else:
        full_mono = mono_ac_a
    if bi_ac is not None:
        mono_and_bi_dict = add_mono_calling_to_hemizygous_loci(bi_ac, full_mono)
        print_ready = flatten_bi_allelic_binning(mono_and_bi_dict)
    else:
        print_ready = full_mono
    ftd = textify_keys_in_mutations_dict(print_ready, sr_label_func, ms_label_func)
    td = filter_bipartition_loci(ftd)
    return td


def filter_chromosomes(mutations_dict, chromosome_biallelic_map):
    fmd = dict()
    by_loc = transpose_dict(mutations_dict)
    for loc in by_loc:
        assert type(loc) == str
        assert loc.split('_')[0] == 'LOC'
        ms = Microsatellite.objects.get(id=int(loc.split('_')[1]))
        if not chromosome_biallelic_map[ms.slice.chromosome.name]:
            continue
        fmd[loc] = by_loc[loc]
    return transpose_dict(fmd)
