import argparse
import clineage.wsgi
import pandas as pd
from sequencing.analysis.full_msv.full_msv import merge
from sequencing.analysis.full_msv.full_msv import split_merged_reads_as_list
from sequencing.analysis.full_msv.full_msv import align_reads_to_ms_variations_part_as_list
from sequencing.analysis.full_msv.full_msv import merge_fmsva_parts
from sequencing.analysis.full_msv.full_msv import separate_reads_by_genotypes
from sequencing.analysis.models import Histogram
from targeted_enrichment.planning.models import Microsatellite
from targeted_enrichment.amplicons.models import Amplicon
from setup.document_targets import document_targets
from setup.document_libraries import document_library
from setup.populate_schemas import a_schema_mono, g_schema_mono, ac_schema_mono, ac_schema_bi
from sequencing.calling.models import BestCorrelationCalledAlleles


def ms_by_amp(ampid):
    l = []
    amp = Amplicon.objects.get(id=ampid)
    for ms in Microsatellite.objects.filter(
        slice__chromosome_id=amp.slice.chromosome_id,
        slice__start_pos__gte=amp.slice.start_pos,
        slice__end_pos__lte=amp.slice.end_pos):
        if ms.microsatellitehistogramgenotype_set.count() == 0:
            continue
        l.append(ms)
    return l


def update_histograms_read_nums():
    for h in Histogram.objects.all():
        if h.num_reads is None:
            num_reads = 0
            for her in h.histogramentryreads_set.all():
                num_reads += her.num_reads
            h.num_reads = num_reads
            h.save()



def collect_histograms_by_type_by_chrom(demultiplexing, min_reads=3):
    ms_by_amp_dict = dict()
    histograms_by_type_by_chrom = dict()
    for h in Histogram.objects.filter(sample_reads__demux=demultiplexing, num_reads__gte=3):
        amp_id = h.amplicon_id
        if amp_id not in ms_by_amp_dict:
            ms_by_amp_dict[amp_id] = ms_by_amp(amp_id)
        for ms in ms_by_amp_dict[amp_id]:
            histograms_by_type_by_chrom.setdefault(ms.repeat_unit_type, dict()).setdefault(ms.slice.chromosome.name,[]).append((h,ms))
    return histograms_by_type_by_chrom


def genotype(histograms_by_type_by_chrom, repeat_unit_type, schema):
    if repeat_unit_type not in histograms_by_type_by_chrom:
        return
    calling_space = [i for l in histograms_by_type_by_chrom[repeat_unit_type].values() for i in l]
    calling_space = [x for x in calling_space if x[0].num_reads >= 3]
    for h, ms in calling_space:
        schema.call_ms_hist(h, ms)


if '__main__' == __name__:
    parser = argparse.ArgumentParser(description='Analyses hist-pairs file')
    parser.add_argument('-i1', '--fastqr1', type=str, dest='fastq_r1', default='input_R1.fastq', help='path of input R1 fastq file')
    parser.add_argument('-i2', '--fastqr2', type=str, dest='fastq_r2', default='input_R2.fastq', help='path of input R2 fastq file')
    parser.add_argument('-r', '--regions', type=str, dest='regions', default='str_regions.bed', help='path of str coordinates bed file')
    parser.add_argument('-o', '--outputfile', type=str, dest='output', default='output.tab', help='path of output file')

    # additional arguments
    parser.add_argument('-mr', '--minreads', type=int, dest='minreads', default=10, help='minimal reads threshold')
    parser.add_argument('-s', '--maxscore', type=float, dest='maxscore', default=999.9, help='max (worst) score allowed')
    parser.add_argument('-v', '--verbose', type=bool, dest='verbose', default=False, help='prints additional calling information columns in the output table')
    parser.add_argument('-ma', '--maxalleles', type=int, dest='max_alleles', default=2, help='maximum number of alleles to consider. *This drastically increase runtime!')
    parser.add_argument('-ml', '--maxmslength', type=int, dest='max_ms_length', default=60, help='maximum ms length to consider')
    parser.add_argument('-md', '--mediandistance', type=int, dest='max_distance_from_median', default=10, help='maximum distance from the median value to take into account')

    args = parser.parse_args()
    fastq_r1 = args.fastq_r1
    fastq_r2 = args.fastq_r2
    regions_path = args.regions
    output_file = args.output

    reads_threshold = args.minreads
    score_threshold = args.maxscore
    max_alleles = args.max_alleles
    max_ms_length = args.max_ms_length
    max_distance_from_median = args.max_distance_from_median
    verbose = args.verbose

    if verbose:
        print('Documenting targets')
    amplicons, ters = document_targets(regions_path)
    if verbose:
        print('Documenting library')
    sr = document_library(amplicons, ters, fastq_r1, fastq_r2)

    if verbose:
        print('Preparing fastqs')
    merged_reads = merge(sr)
    fmsv_merged_reads_parts_lists = split_merged_reads_as_list(merged_reads,
                                                               reads_chunk_size=10 ** 5, included_reads='M'
                                                               )
    if verbose:
        print('Aligning reads')
    fmsva_parts_lists = align_reads_to_ms_variations_part_as_list(
        fmsv_merged_reads_parts_lists[0], padding=50, mss_version=0, chunk_size=15000)
    merged_fmsvas = merge_fmsva_parts(
        [fmsva_parts_lists],
        reads_chunk_size=10 ** 5,
        included_reads='M'
    )
    l = list(merged_fmsvas)
    if verbose:
        print('Generating Histograms')
    fhers_list = separate_reads_by_genotypes(l[0])
    l2 = list(fhers_list)
    update_histograms_read_nums()
    if verbose:
        print('Collecting Histograms')
    histograms_by_type_by_chrom = collect_histograms_by_type_by_chrom(sr.demux, min_reads=reads_threshold)

    if verbose:
        print('Genotyping')
    genotype(histograms_by_type_by_chrom, 'A', a_schema_mono)
    genotype(histograms_by_type_by_chrom, 'C', g_schema_mono)
    genotype(histograms_by_type_by_chrom, 'AC', ac_schema_mono)
    try:
        genotype(histograms_by_type_by_chrom, 'AC', ac_schema_bi)
    except AssertionError:
        print('Problem with bi-allelic genotyping')

    if verbose:
        print('Formatting')
    d = dict()
    for ca in BestCorrelationCalledAlleles.objects.filter(
        histogram__in=Histogram.objects.filter(
            num_reads__gte=reads_threshold,
            sample_reads=sr),
        confidence__lte=score_threshold).select_related(
            'histogram__sample_reads', 'microsatellite__slice__chromosome', 'genotypes'):
        d[ca.microsatellite.id] = {
            'Name': ca.microsatellite.name,
            'Genotype': ca.genotypes.alleles,
            'Confidence': ca.confidence,
            'Cycle': ca.cycle,
            'Scheme': ca.calling_scheme.name,
        }

    df = pd.DataFrame.from_dict(d, orient="index")
    df.to_csv(output_file)

