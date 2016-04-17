
import itertools
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
    ms_genotypes_to_name, AdamHistogram, HistogramEntryReads
from targeted_enrichment.planning.models import Microsatellite

pear = local["pear"]
pear_with_defaults = pear["-v", "40",
                          "-m", "300"]


bowtie2build = local["bowtie2-build"]

# bowtie2-build \
#     ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged.fa \
#     ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged

bowtie2 = local["bowtie2"]
bowtie2_with_defaults = bowtie2["-p", "24",
                                "-a",
                                "--very-sensitive",]
bowtie2_with_defaults2 = bowtie2["-p", "24",
                                 "-a"]

# bowtie2 \
#     -p 24 \
#     -a \
#     --very-sensitive \
#     -x ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged \
#     -f ${TargetsFolder}/${TargetsName}_Primers.fa \
#         > ${fastq_folder}/Cluster/processing/${fastqName}/merged/${fastq_R1}_final_merged_Primers.sam

# /net/mraid11/export/data/dcsoft/home/Adam/Software/bowtie2-2.2.2/bowtie2 \
#     -a --very-sensitive \
#     -x index_1 \
#     -f /net/mraid11/export/dcstor/Ofir/ngs_fixtures/28727_and_28734_Primers.fa \
#     -S test2.sam



#~/data/home/Adam/Software/PEAR-master/src/pear
    # -f  ${fastq_folder}/Cluster/processing/${fastqName}/${fastq_R1}_cutadapt_ListBoth.fastq
    # -r ${fastq_folder}/Cluster/processing/${fastqName}/${fastq_R2}_cutadapt_ListBoth.fastq
    # -o ${fastq_folder}/Cluster/processing/${fastqName}/${fastq_R1}_cutadapt_ListBoth_PEAR
    # -v 40
    # -m 300


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
    fasta = _create_padded_fasta(reads, padding)
    # TODO: clean this double code.
    os.mkdir(index_dir)
    ri = AdamReadsIndex(
        merged_reads=merged_reads,
        included_reads=included_reads,
        index_dump_dir=index_dir,
        padding=padding,
    )
    bowtie2build(fasta, ri.index_files_prefix)
    os.unlink(fasta)
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
    amplicons = reads_index.merged_reads.sample_reads.library.amplicons
    panel_fasta = _create_panel_fasta(amplicons)
    bowtie2_with_defaults('-x', reads_index.index_files_prefix,
                          '-f', panel_fasta,
                          '-S', assignment_sam)
    os.unlink(panel_fasta)
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


def _extract_amplicon_reads(indexed_reads, read_ids):
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
        amplicon_readsm_fastq_name = _extract_amplicon_reads(reads, read_ids)
        amplicon_reads1_fastq_name = _extract_amplicon_reads(reads1, read_ids)
        amplicon_reads2_fastq_name = _extract_amplicon_reads(reads2, read_ids)
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


def _build_ms_variations(amplicon, padding):
    amplicon = amplicon.subclass
    # FIXME: What's the right API for this?
    targets = list(amplicon.ter.te.targets.select_related(
        "microsatellite",
        "slice",
    ))
    mss = []
    for t in targets:
        try:
            m = t.microsatellite
        except Microsatellite.DoesNotExist:
            pass
        else:
            mss.append(m)
    # This is so they are ordered properly for the format string.
    mss = sorted(mss,key=lambda ms: ms.slice.start_pos)
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


def get_adam_ms_variations(amplicon, padding):
    try:
        return AdamMSVariations.objects.get(
            amplicon=amplicon,
            padding=padding,
        )
    except AdamMSVariations.DoesNotExist:
        index_dir = get_unique_path()
        fasta = _build_ms_variations(amplicon, padding)
        os.mkdir(index_dir)
        mock_msv = AdamMSVariations(
            amplicon=amplicon,
            padding=padding,
            index_dump_dir=index_dir,
        )
        bowtie2build(fasta, mock_msv.index_files_prefix)
        os.unlink(fasta)
        msv, c = AdamMSVariations.objects.get_or_create(
            amplicon=amplicon,
            padding=padding,
            defaults={"index_dump_dir": index_dir}
        )
        if not c:
            delete_files(mock_msv)
        return msv


def align_reads_to_ms_variations(amplicon_reads, padding):
    assignment_sam = get_unique_path("sam")
    msv = get_adam_ms_variations(amplicon_reads.amplicon, padding)
    bowtie2_with_defaults2('-x', msv.index_files_prefix,
                          '-U', amplicon_reads.fastqm,
                          '-S', assignment_sam)
    ah = AdamHistogram.objects.create(
        sample_reads_id=amplicon_reads.margin_assignment.reads_index \
            .merged_reads.sample_reads_id,
        amplicon_reads=amplicon_reads,
        assignment_sam=assignment_sam,
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
        genotypes_readsm_fastq_name = _extract_amplicon_reads(readsm, read_ids)
        genotypes_reads1_fastq_name = _extract_amplicon_reads(reads1, read_ids)
        genotypes_reads2_fastq_name = _extract_amplicon_reads(reads2, read_ids)
        her = HistogramEntryReads.objects.create(
            histogram=histogram,
            amplicon=histogram.amplicon_reads.amplicon,
            num_reads=len(read_ids),
            fastq1=genotypes_reads1_fastq_name,
            fastq2=genotypes_reads2_fastq_name,
            fastqm=genotypes_readsm_fastq_name,
        )
        her.microsatellite_genotypes.add(*genotypes)
        yield her
