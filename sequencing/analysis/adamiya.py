
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

from misc.utils import get_unique_path
from sequencing.analysis.models import AdamMergedReads, AdamReadsIndex, \
    AdamMarginAssignment, AdamAmpliconReads, amplicon_margin_to_name, \
    LEFT, RIGHT, AdamMSVariations, MicrosatelliteHistogramGenotype, \
    ms_genotypes_to_name, AdamHistogram, HistogramEntryReads, SampleReads, \
    delete_files
from targeted_enrichment.planning.models import Microsatellite

from distributed.executor import as_completed

pear = local["pear"]
pear_with_defaults = pear["-v", "40",
                          "-m", "300"]


bowtie2build = local["bowtie2-build"]
bowtie2build_fixed_seed = bowtie2build["--seed", "1"]

bowtie2 = local["bowtie2"]
bowtie2_fixed_seed = bowtie2["--seed", "1"]
bowtie2_with_defaults = bowtie2_fixed_seed["-p", "24",
                                "-a",
                                "--very-sensitive"]
bowtie2_with_defaults2 = bowtie2_fixed_seed["-p", "24",
                                 "-a"]


@contextlib.contextmanager
def unlink(fname):
    try:
        yield fname
    finally:
        os.unlink(fname)


def merge(sample_reads):
    prefix = get_unique_path()
    # TODO: check by merging scheme.
    pear_with_defaults("-f", sample_reads.fastq1,
                       "-r", sample_reads.fastq2,
                       "-o", prefix)
    mr = AdamMergedReads.objects.create(
        sample_reads=sample_reads,
        assembled_fastq="{}.assembled.fastq".format(prefix),
        discarded_fastq="{}.discarded.fastq".format(prefix),
        unassembled_forward_fastq="{}.unassembled.forward.fastq".format(prefix),
        unassembled_reverse_fastq="{}.unassembled.reverse.fastq".format(prefix),
    )
    return mr

def _pad_records(records, padding):
    for record in records:
        yield "N"*padding + record + "N"*padding

def _create_padded_fasta(reads_iter, padding):
    file_path = get_unique_path("fasta")
    padded_reads = _pad_records(reads_iter, padding)
    SeqIO.write(padded_reads, file_path, "fasta")
    return file_path

def create_reads_index(merged_reads, included_reads, padding):
    """
    Create an index from some of these reads. included_reads should be one
    of ReadsIndex.INCLUDED_READS_OPTIONS, and chooses which of the reads we take.
    padding controls how much to pad on each side of the reads.
    """
    index_dir = get_unique_path()
    reads = merged_reads.included_reads_generator(included_reads)
    with unlink(_create_padded_fasta(reads, padding)) as fasta:
        # TODO: clean this double code.
        os.mkdir(index_dir)
        ri = AdamReadsIndex(
            merged_reads=merged_reads,
            included_reads=included_reads,
            index_dump_dir=index_dir,
            padding=padding,
        )
        bowtie2build_fixed_seed(fasta, ri.index_files_prefix)
    ri.save()
    return ri


def _primers_seqrecords_generator(amplicons):
    for uw in amplicons:
        # NOTE: [1:] is workaround adam bug ?
        yield SeqRecord(Seq(str(uw.left_margin)), id=amplicon_margin_to_name(uw, LEFT), name='', description='')[1:]
        yield SeqRecord(Seq(str(uw.right_margin)), id=amplicon_margin_to_name(uw, RIGHT), name='', description='')[1:]
        # name='', description='' are workarounds for the '<unknown description>' that is being outputted and
        # breaks the downstream bowtie2 alignment.


def _create_panel_fasta(amplicons):
    panel_fasta_name = get_unique_path("fasta")
    SeqIO.write(_primers_seqrecords_generator(amplicons),
        panel_fasta_name, "fasta")
    return panel_fasta_name


def align_primers_to_reads(reads_index):
    assignment_sam = get_unique_path("sam")
    amplicons = reads_index.merged_reads.sample_reads.library.subclass.amplicons
    with unlink(_create_panel_fasta(amplicons)) as panel_fasta:
        bowtie2_with_defaults('-x', reads_index.index_files_prefix,
                              '-f', panel_fasta,
                              '-S', assignment_sam)
    ama = AdamMarginAssignment.objects.create(
        reads_index=reads_index,
        assignment_sam=assignment_sam,
    )
    return ama


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


def _extract_reads_by_id(indexed_reads, read_ids):
    amplicon_reads_fastq_name = get_unique_path("fastq")
    amplicon_reads = (indexed_reads[read_id] for read_id in read_ids)
    SeqIO.write(amplicon_reads, amplicon_reads_fastq_name, "fastq")
    return amplicon_reads_fastq_name


def seperate_reads_by_amplicons(margin_assignment):
    reads_matches = _collect_mappings_from_sam(margin_assignment)
    validated_reads_amplicons = _validate_amplicon_mapping(reads_matches)
    reads_by_amplicon = _aggregate_read_ids_by_amplicon(validated_reads_amplicons)
    reads_gen = margin_assignment.reads_index.included_reads_generator()
    reads1 = SeqIO.index(margin_assignment.reads_index.merged_reads \
        .sample_reads.fastq1, "fastq")
    reads2 = SeqIO.index(margin_assignment.reads_index.merged_reads \
        .sample_reads.fastq2, "fastq")
    reads = SeqIO.to_dict(reads_gen)
    for amplicon, read_ids in reads_by_amplicon.items():
        amplicon_readsm_fastq_name = _extract_reads_by_id(reads, read_ids)
        amplicon_reads1_fastq_name = _extract_reads_by_id(reads1, read_ids)
        amplicon_reads2_fastq_name = _extract_reads_by_id(reads2, read_ids)
        aar = AdamAmpliconReads.objects.create(
            margin_assignment=margin_assignment,
            amplicon=amplicon,
            fastqm=amplicon_readsm_fastq_name,
            fastq1=amplicon_reads1_fastq_name,
            fastq2=amplicon_reads2_fastq_name,
        )
        yield aar


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
    fasta = get_unique_path("fa")
    prefix = "{}".format(amplicon.id)
    SeqIO.write(_get_mss_variations_seqrecords(mss, full_fmt, prefix),
        fasta, "fasta")
    return fasta


def get_adam_ms_variations(amplicon, padding, mss_version):
    try:
        return AdamMSVariations.objects.get(
            amplicon=amplicon,
            padding=padding,
            microsatellites_version=mss_version,
        )
    except AdamMSVariations.DoesNotExist:
        index_dir = get_unique_path()
        mss = Microsatellite.objects.filter(
            slice__start_pos__gte=amplicon.slice.start_pos,
            slice__end_pos__lte=amplicon.slice.end_pos,
            slice__chromosome_id=amplicon.slice.chromosome_id,
            planning_version=mss_version,
        )
        with unlink(_build_ms_variations(amplicon, padding, mss)) as fasta:
            os.mkdir(index_dir)
            mock_msv = AdamMSVariations(
                amplicon=amplicon,
                padding=padding,
                microsatellites_version=mss_version,
                index_dump_dir=index_dir,
            )
            bowtie2build(fasta, mock_msv.index_files_prefix)
        msv, c = AdamMSVariations.objects.get_or_create(
            amplicon=amplicon,
            padding=padding,
            microsatellites_version=mss_version,
            defaults={"index_dump_dir": index_dir}
        )
        if not c:
            delete_files(AdamMSVariations, mock_msv)
        return msv


def align_reads_to_ms_variations(amplicon_reads, padding, mss_version):
    assignment_sam = get_unique_path("sam")
    msv = get_adam_ms_variations(amplicon_reads.amplicon.subclass, padding, 
        mss_version)
    bowtie2_with_defaults2('-x', msv.index_files_prefix,
                          '-U', amplicon_reads.fastqm,
                          '-S', assignment_sam)
    ah = AdamHistogram.objects.create(
        sample_reads_id=amplicon_reads.margin_assignment.reads_index \
            .merged_reads.sample_reads_id,
        amplicon_reads=amplicon_reads,
        assignment_sam=assignment_sam,
        ms_variations=msv,
    )
    return ah


def _collect_genotypes_from_sam(histogram):
    genotypes_reads = set()
    for read_id, ms_genotypes in histogram.read_sam():
        if read_id not in genotypes_reads:
            # We assume that the SAM is ordered by quality of match.
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
    for genotypes, read_ids in genotypes_reads.items():
        genotypes_readsm_fastq_name = _extract_reads_by_id(readsm, read_ids)
        genotypes_reads1_fastq_name = _extract_reads_by_id(reads1, read_ids)
        genotypes_reads2_fastq_name = _extract_reads_by_id(reads2, read_ids)
        her = HistogramEntryReads.objects.create(
            histogram=histogram,
            amplicon=histogram.amplicon_reads.amplicon,
            microsatellites_version=histogram.ms_variations \
                .microsatellites_version,
            num_reads=len(read_ids),
            fastq1=genotypes_reads1_fastq_name,
            fastq2=genotypes_reads2_fastq_name,
            fastqm=genotypes_readsm_fastq_name,
        )
        her.microsatellite_genotypes.add(*genotypes)
        yield her


def double_map(executor, func, future_lists, *params):
    for f in as_completed(future_lists):
        l = f.result()
        yield executor.map(func, l, *[itertools.repeat(p) for p in params], pure=False)


def list_iterator(f):
    @functools.wraps(f)
    def inner(*args, **kwargs):
        return list(f(*args, **kwargs))
    return inner


def run_parallel(executor, sample_reads, included_reads="F", mss_version=0, read_padding=5, ref_padding=50):
    merged_reads = executor.map(merge, sample_reads, pure=False)
    reads_indices = executor.map(create_reads_index, merged_reads,
        itertools.repeat(included_reads), itertools.repeat(read_padding), pure=False)
    adam_margin_assignments = executor.map(align_primers_to_reads,
        reads_indices, pure=False)
    adam_amplicon_reads_lists = executor.map(
        list_iterator(seperate_reads_by_amplicons), adam_margin_assignments, pure=False)
    yield merged_reads, reads_indices, adam_margin_assignments, adam_amplicon_reads_lists
    adam_histograms = double_map(executor, align_reads_to_ms_variations,
        adam_amplicon_reads_lists, ref_padding, mss_version)
    for fs in adam_histograms:
        fs2 = executor.map(list_iterator(separate_reads_by_genotypes), fs, pure=False)
        yield fs, fs2
