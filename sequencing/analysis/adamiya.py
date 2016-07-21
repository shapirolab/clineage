
import itertools
import functools
import contextlib
from plumbum import local
import os
import uuid
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil

from django.conf import settings
from django.db import IntegrityError

from misc.utils import get_unique_path, unique_file_cm, unique_dir_cm, unlink, \
    raise_or_create, get_get_or_create
from sequencing.analysis.models import AdamMergedReads, AdamReadsIndex, \
    AdamMarginAssignment, AdamAmpliconReads, amplicon_margin_to_name, \
    LEFT, RIGHT, AdamMSVariations, MicrosatelliteHistogramGenotype, \
    MicrosatelliteHistogramGenotypeSet, ms_genotypes_to_name, AdamHistogram, \
    HistogramEntryReads, SampleReads, delete_files, name_to_ms_genotypes, \
    BowtieIndexMixin, PearOutputMixin, SNPHistogramGenotypeSet, SNPHistogramGenotype
from targeted_enrichment.planning.models import Microsatellite

from distributed.executor import as_completed

pear = local["pear"]
pear_with_defaults = pear["-v", "40",
                          "-m", "300"]


bowtie2build = local["bowtie2-build"]
bowtie2build_fixed_seed = bowtie2build["--seed", "1"]

bowtie2 = local["bowtie2"]
bowtie2_fixed_seed = bowtie2["--seed", "1"]
bowtie2_with_defaults = bowtie2_fixed_seed["-p", "1",
                                "-a",
                                "--very-sensitive"]
bowtie2_with_defaults2 = bowtie2_fixed_seed["-p", "1",
                                 "-a"]


def merge(sample_reads):
    def inner():
        with unique_dir_cm() as pear_dir:
            pear_output = PearOutputMixin(
                pear_dump_dir=pear_dir,
            )
            pear_with_defaults("-f", sample_reads.fastq1,
                               "-r", sample_reads.fastq2,
                               "-o", pear_output.pear_files_prefix)
            return raise_or_create(AdamMergedReads,
                sample_reads=sample_reads,
                defaults=dict(
                    pear_dump_dir=pear_dir,
                ),
            )
    return get_get_or_create(inner, AdamMergedReads,
        sample_reads=sample_reads,
    )


def _pad_records(records, padding):
    for record in records:
        yield "N"*padding + record + "N"*padding


def _create_padded_fasta(reads_iter, padding):
    padded_reads = _pad_records(reads_iter, padding)
    with unique_file_cm("fasta") as file_path:
        SeqIO.write(padded_reads, file_path, "fasta")
        return file_path


def create_reads_index(merged_reads, included_reads, padding):
    """
    Create an index from some of these reads. included_reads should be one
    of ReadsIndex.INCLUDED_READS_OPTIONS, and chooses which of the reads we take.
    padding controls how much to pad on each side of the reads.
    """
    def inner():
        reads = merged_reads.included_reads_generator(included_reads)
        with unique_dir_cm() as index_dir:
            bowtie_index = BowtieIndexMixin(index_dump_dir=index_dir)
            with unlink(_create_padded_fasta(reads, padding)) as fasta:
                bowtie2build_fixed_seed(fasta, bowtie_index.index_files_prefix)
            return raise_or_create(AdamReadsIndex,
                merged_reads=merged_reads,
                included_reads=included_reads,
                padding=padding,
                defaults=dict(
                    index_dump_dir=index_dir,
                ),
            )
    return get_get_or_create(inner, AdamReadsIndex,
        merged_reads=merged_reads,
        included_reads=included_reads,
        padding=padding,
    )


def _primers_seqrecords_generator(amplicons):
    for uw in amplicons:
        # NOTE: [1:] is workaround adam bug ?
        yield SeqRecord(Seq(str(uw.left_margin)), id=amplicon_margin_to_name(uw, LEFT), name='', description='')[1:]
        yield SeqRecord(Seq(str(uw.right_margin)), id=amplicon_margin_to_name(uw, RIGHT), name='', description='')[1:]
        # name='', description='' are workarounds for the '<unknown description>' that is being outputted and
        # breaks the downstream bowtie2 alignment.


def _create_panel_fasta(amplicons):
    panel_primers = _primers_seqrecords_generator(amplicons)
    with unique_file_cm("fasta") as panel_fasta_name:
        SeqIO.write(panel_primers, panel_fasta_name, "fasta")
        return panel_fasta_name


def align_primers_to_reads(reads_index):
    def inner():
        amplicons = reads_index.merged_reads.sample_reads.library \
            .subclass.amplicons
        with unique_file_cm("sam") as assignment_sam:
            with unlink(_create_panel_fasta(amplicons)) as panel_fasta:
                bowtie2_with_defaults('-x', reads_index.index_files_prefix,
                                      '-f', panel_fasta,
                                      '-S', assignment_sam)
            return raise_or_create(AdamMarginAssignment,
                reads_index=reads_index,
                defaults=dict(
                    assignment_sam=assignment_sam,
                ),
            )
    return get_get_or_create(inner, AdamMarginAssignment,
        reads_index=reads_index)


def _collect_mappings_from_sam(margin_assignment):
    reads_matches = {}
    for read_id, amplicon_margin in margin_assignment.read_sam():
        reads_matches.setdefault(read_id, set()).add(amplicon_margin)
    return reads_matches


def _validate_amplicon_mapping(reads_matches):
    # TODO: consider maintaining pysqm alignment objects for stronger tests here.
    # For example, make sure left primer is mapped to the left of the amplicon
    # and the right primer to the right of the amplicon using ".pos" property
    # of pysam.calignmentfile.AlignedSegment objects
    for read_id, matches in reads_matches.items():
        if len(matches) != 2:
            # place proper bin
            continue
        amplicons, directions = list(zip(*matches))
        u1, u2 = amplicons
        if u1 != u2:
            # place proper bin
            continue
        if set(directions) != {LEFT, RIGHT}:
            raise RuntimeError("What else is here? {}".format(directions))
        yield read_id, u1


def _aggregate_read_ids_by_amplicon(validated_reads_amplicons):
    reads_by_amplicon = {}
    for read_id, amplicon in validated_reads_amplicons:
        reads_by_amplicon.setdefault(amplicon, set()).add(read_id)
    return reads_by_amplicon


@contextlib.contextmanager
def _extract_reads_by_id(indexed_reads, read_ids):
    amplicon_reads = (indexed_reads[read_id] for read_id in read_ids)
    with unique_file_cm("fastq") as amplicon_reads_fastq_name:
        SeqIO.write(amplicon_reads, amplicon_reads_fastq_name, "fastq")
        yield amplicon_reads_fastq_name


def separate_reads_by_amplicons(margin_assignment):
    reads_matches = _collect_mappings_from_sam(margin_assignment)
    validated_reads_amplicons = _validate_amplicon_mapping(reads_matches)
    reads_by_amplicon = _aggregate_read_ids_by_amplicon( \
        validated_reads_amplicons)
    reads_gen = margin_assignment.reads_index.included_reads_generator()
    reads1 = SeqIO.index(margin_assignment.reads_index.merged_reads \
        .sample_reads.fastq1, "fastq")
    reads2 = SeqIO.index(margin_assignment.reads_index.merged_reads \
        .sample_reads.fastq2, "fastq")
    reads = SeqIO.to_dict(reads_gen)
    for amplicon, read_ids in reads_by_amplicon.items():
        def inner():
            with _extract_reads_by_id(reads, read_ids) as \
                    amplicon_readsm_fastq_name, \
                _extract_reads_by_id(reads1, read_ids) as \
                    amplicon_reads1_fastq_name, \
                _extract_reads_by_id(reads2, read_ids) as \
                    amplicon_reads2_fastq_name:
                return raise_or_create(AdamAmpliconReads,
                    margin_assignment=margin_assignment,
                    amplicon=amplicon,
                    defaults=dict(
                        fastqm=amplicon_readsm_fastq_name,
                        fastq1=amplicon_reads1_fastq_name,
                        fastq2=amplicon_reads2_fastq_name,
                    ),
                )
        yield get_get_or_create(inner, AdamAmpliconReads, 
            margin_assignment=margin_assignment,
            amplicon=amplicon,
        )


def _get_ms_length_range(ms):
    #if (Rep>3)
            #Xcomb{ti} = 3:(Rep*2-3);
    #else
            #Xcomb{ti} = 1:5;
    #end
    # TODO: put better boundries
    n = int(ms.repeat_number)
    if n > 3:
        return range(3,2*n-2)
    else:
        return range(1,6)


def _get_ms_variations(ms):
    for i in _get_ms_length_range(ms):
        yield MicrosatelliteHistogramGenotype.get_for_genotype(ms,i)


def _get_mss_variations_seqrecords(mss, seq_fmt, prefix):
    for mult in itertools.product(*[_get_ms_variations(ms) for ms in mss]):
        seqs = [mhg.sequence for mhg in mult]
        # name='', description='' are workarounds for the '<unknown
        # description>' that is being outputted otherwise
        yield SeqRecord(
            Seq(seq_fmt.format(*seqs)),
            id=ms_genotypes_to_name(mult, prefix),
            name='',
            description=''
        )


def _build_ms_variations(amplicon, padding, mss):
    amplicon = amplicon.subclass
    # This is so they are ordered properly for the format string.
    mss = sorted(mss, key=lambda ms: ms.slice.start_pos)
    # FIXME: kill this +-1 when we move to 0-based.
    points = [amplicon.slice.start_pos-1]
    for ms in mss:
        # FIXME: kill this +-1 when we move to 0-based.
        points.append(ms.slice.start_pos-1)
        points.append(ms.slice.end_pos)
    points.append(amplicon.slice.end_pos)
    if points != sorted(points):
        raise IntegrityError("Amplicon {} has interlocking MSs or MSs " \
            "outside its boundaries".format(amplicon.id))
    # FIXME: kill this +-1 when we move to 0-based.
    segments = [(points[2*i]+1, points[2*i+1]) \
        for i in range(len(points)//2) if points[2*i] < points[2*i+1]]
    # FIXME: maybe store these slices?
    fmt = "{}".join([amplicon.slice.chromosome.getdna(*a) for a in segments])
    full_fmt = "{pad}{left}{slice}{right}{pad}".format(
        pad="N"*padding,
        left=amplicon.left_margin,
        slice=fmt,
        right=amplicon.right_margin
    )
    prefix = "{}".format(amplicon.id)
    with unique_file_cm("fa") as fasta:
        SeqIO.write(_get_mss_variations_seqrecords(mss, full_fmt, prefix),
            fasta, "fasta")
        return fasta


def get_adam_ms_variations(amplicon, padding, mss_version):
    def inner():
        mss = Microsatellite.objects.filter(
            slice__start_pos__gte=amplicon.slice.start_pos,
            slice__end_pos__lte=amplicon.slice.end_pos,
            slice__chromosome_id=amplicon.slice.chromosome_id,
            planning_version=mss_version,
        )
        with unique_dir_cm() as index_dir:
            bowtie_index = BowtieIndexMixin(index_dump_dir=index_dir)
            with unlink(_build_ms_variations(amplicon, padding, mss)) as fasta:
                bowtie2build(fasta, bowtie_index.index_files_prefix)
            return raise_or_create(AdamMSVariations,
                amplicon=amplicon,
                padding=padding,
                microsatellites_version=mss_version,
                defaults=dict(
                    index_dump_dir=index_dir,
                ),
            )
    return get_get_or_create(inner, AdamMSVariations,
        amplicon=amplicon,
        padding=padding,
        microsatellites_version=mss_version,
    )


def align_reads_to_ms_variations(amplicon_reads, padding, mss_version):
    msv = get_adam_ms_variations(amplicon_reads.amplicon.subclass, padding,
        mss_version)
    def inner():
        with unique_file_cm("sam") as assignment_sam:
            bowtie2_with_defaults2('-x', msv.index_files_prefix,
                                   '-U', amplicon_reads.fastqm,
                                   '-S', assignment_sam)
            return raise_or_create(AdamHistogram,
                amplicon_reads=amplicon_reads,
                ms_variations=msv,
                defaults=dict(
                    assignment_sam=assignment_sam,
                    sample_reads_id=amplicon_reads.margin_assignment \
                        .reads_index.merged_reads.sample_reads_id,
                    amplicon=amplicon_reads.amplicon,
                    microsatellites_version=msv.microsatellites_version,
                ),
            )
    return get_get_or_create(inner, AdamHistogram,
        amplicon_reads=amplicon_reads,
        ms_variations=msv,
    )


def _collect_genotypes_from_sam(histogram):
    genotypes_reads = set()
    for read_id, ms_genotypes_name in histogram.read_sam():
        if read_id not in genotypes_reads:
            # We assume that the SAM is ordered by quality of match.
            ms_genotypes, prefix = name_to_ms_genotypes(ms_genotypes_name)
            assert int(prefix) == histogram.amplicon_reads.amplicon_id
            genotypes_reads.add(read_id)
            yield ms_genotypes, read_id


def _aggregate_read_ids_by_genotypes(genotype_read_iterator):
    genotypes_reads = {}
    for ms_genotypes, read_id in genotype_read_iterator:
        genotypes_reads.setdefault(ms_genotypes, set()).add(read_id)
    return genotypes_reads


def separate_reads_by_genotypes(histogram):
    genotype_read_iterator = _collect_genotypes_from_sam(histogram)
    genotypes_reads = _aggregate_read_ids_by_genotypes(genotype_read_iterator)
    reads1 = SeqIO.index(histogram.amplicon_reads.fastq1, "fastq")
    reads2 = SeqIO.index(histogram.amplicon_reads.fastq2, "fastq")
    readsm = SeqIO.index(histogram.amplicon_reads.fastqm, "fastq")
    none_ms_genotype = MicrosatelliteHistogramGenotype.objects.get(microsatellite=None)
    none_snp_genotype = SNPHistogramGenotype.objects.get(snp=None)
    snp_histogram_genotypes, c = SNPHistogramGenotypeSet.objects.get_or_create(
        **{fn: none_snp_genotype for fn in SNPHistogramGenotypeSet.genotype_field_names()})
    for genotypes, read_ids in genotypes_reads.items():
        ordered_genotypes = dict(itertools.zip_longest(
            MicrosatelliteHistogramGenotypeSet.genotype_field_names(),
            sorted(list(genotypes), key=lambda g: g.microsatellite.slice),
            fillvalue=none_ms_genotype,
        ))
        microsatellite_histogram_genotypes, c = MicrosatelliteHistogramGenotypeSet.objects.get_or_create(**ordered_genotypes)
        def inner():
            with _extract_reads_by_id(readsm, read_ids) as genotypes_readsm_fastq_name, \
                _extract_reads_by_id(reads1, read_ids) as genotypes_reads1_fastq_name, \
                _extract_reads_by_id(reads2, read_ids) as genotypes_reads2_fastq_name:
                return raise_or_create(HistogramEntryReads,
                    histogram=histogram,
                    microsatellite_genotypes=microsatellite_histogram_genotypes,
                    snp_genotypes=snp_histogram_genotypes,
                    defaults=dict(
                        num_reads=len(read_ids),
                        fastq1=genotypes_reads1_fastq_name,
                        fastq2=genotypes_reads2_fastq_name,
                        fastqm=genotypes_readsm_fastq_name,
                    ),
                )
        yield get_get_or_create(inner, HistogramEntryReads, 
            histogram=histogram,
            microsatellite_genotypes=microsatellite_histogram_genotypes,
            snp_genotypes=snp_histogram_genotypes,
        )


def close_connection_and(f):
    @functools.wraps(f)
    def inner(*args, **kwargs):
        from django.db import connection
        connection.close()
        return f(*args, **kwargs)
    return inner


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
