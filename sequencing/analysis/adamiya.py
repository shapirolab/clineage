
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
    AdamMarginAssignment, AdamAmpliconReads, unwrapper_margin_to_name, \
    LEFT, RIGHT, AdamMSVariations, MicrosatelliteHistogramGenotype
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
        demux_reads=sample_reads,
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


def _primers_seqrecords_generator(unwrappers):
    for uw in unwrappers:
        # NOTE: [1:] is workaround adam bug ?
        yield SeqRecord(Seq(uw.left_margin.seq), id=unwrapper_margin_to_name(uw, LEFT), name='', description='')[1:]
        yield SeqRecord(Seq(uw.right_margin.seq), id=unwrapper_margin_to_name(uw, RIGHT), name='', description='')[1:]
        # name='', description='' are workarounds for the '<unknown description>' that is being outputted and
        # breaks the downstream bowtie2 alignment.


def _create_panel_fasta(unwrappers):
    panel_fasta_name = get_unique_path("fasta")
    SeqIO.write(_primers_seqrecords_generator(unwrappers),
        panel_fasta_name, "fasta")
    return panel_fasta_name


def align_primers_to_reads(reads_index):
    assignment_sam = get_unique_path("sam")
    unwrappers = reads_index.merged_reads.demux_reads.library.unwrappers
    panel_fasta = _create_panel_fasta(unwrappers)
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
    for read_id, unwrapper_margin in margin_assignment.read_sam():
        reads_matches.setdefault(read_id, set()).add(unwrapper_margin)
    return reads_matches


def _validate_unwrapper_mapping(reads_matches):
    # TODO: consider maintaining pysqm alignment objects for stronger tests here.
    # For example, make sure left primer is mapped to the left of the amplicon
    # and the right primer to the right of the amplicon using ".pos" property
    # of pysam.calignmentfile.AlignedSegment objects
    for read_id, matches in reads_matches.iteritems():
        if len(matches) != 2:
            # place proper bin
            continue
        unwrappers, directions = zip(*matches)
        u1, u2 = unwrappers
        if u1 != u2:
            # place proper bin
            continue
        if set(directions) != {LEFT, RIGHT}:
            raise RuntimeError("What else is here? {}".format(directions))
        yield read_id, u1


def _aggregate_read_ids_by_unwrapper(validated_reads_unwrappers):
    reads_by_unwrapper = {}
    for read_id, unwrapper in validated_reads_unwrappers:
        reads_by_unwrapper.setdefault(unwrapper, set()).add(read_id)
    return reads_by_unwrapper


def seperate_reads_by_amplicons(margin_assignment):
    reads_matches = _collect_mappings_from_sam(margin_assignment)
    validated_reads_unwrappers = _validate_unwrapper_mapping(reads_matches)
    reads_by_unwrapper = _aggregate_read_ids_by_unwrapper(validated_reads_unwrappers)
    reads = margin_assignment.reads_index.included_reads_generator()
    fq_dict = SeqIO.to_dict(reads)
    for unwrapper, read_ids in reads_by_unwrapper.iteritems():
        unwrapper_reads_fastq_name = get_unique_path("fastq")
        unwrapper_reads = (fq_dict[read_id] for read_id in read_ids)
        SeqIO.write(unwrapper_reads, unwrapper_reads_fastq_name, "fastq")
        aar = AdamAmpliconReads.objects.create(
            margin_assignment=margin_assignment,
            unwrapper=unwrapper,
            fastq=unwrapper_reads_fastq_name,
        )
        yield aar


def _get_ms_length_range(ms):
    #if (Rep>3)
            #Xcomb{ti} = 3:(Rep*2-3);
    #else
            #Xcomb{ti} = 1:5;
    #end
    # TODO: put better boundries
    n = ms.repeat_number
    if n > 3:
        return xrange(3,2*n-2)
    else:
        return xrange(1,6)


def _get_ms_variations(ms):
    for i in _get_ms_length_range(ms):
        yield MicrosatelliteHistogramGenotype.get_for_genotype(ms,i)


def _get_mss_variations_seqrecords(mss, seq_fmt, prefix):
    for mult in itertools.product(*[_get_ms_variations(ms) for ms in mss]):
        seqs = [mhg.sequence for mhg in mult]
        names = ["{}".format(mhg) for mhg in mult]
        # name='', description='' are workarounds for the '<unknown
        # description>' that is being outputted otherwise
        yield SeqRecord(
            Seq(seq_fmt.format(*seqs)),
            id=":".join([prefix] + names),
            name='',
            description=''
        )


def _build_ms_variations(unwrapper, padding):
    # FIXME: What's the right API for this?
    targets = list(unwrapper.ter.te.targets.select_related(
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
    points = [unwrapper.slice.start_pos-1]
    for ms in mss:
        # FIXME: kill this +-1 when we move to 0-based.
        points.append(ms.slice.start_pos-1)
        points.append(ms.slice.end_pos)
    points.append(unwrapper.slice.end_pos)
    if points != sorted(points):
        raise IntegrityError("Unwrapper {} has interlocking MSs or MSs " \
            "outside its boundaries".format(unwrapper.id))
    # FIXME: kill this +-1 when we move to 0-based.
    segments = [(points[2*i]+1, points[2*i+1]) \
        for i in xrange(len(points)//2) if points[2*i] < points[2*i+1]]
    # FIXME: maybe store these slices?
    fmt = "{}".join([unwrapper.slice.chromosome.getdna(*a) for a in segments])
    full_fmt = "{pad}{left}{slice}{right}{pad}".format(
        pad="N"*padding,
        left=unwrapper.left_margin,
        slice=fmt,
        right=unwrapper.right_margin
    )
    fasta = get_unique_path("fa")
    prefix = "{}".format(unwrapper.id)
    SeqIO.write(_get_mss_variations_seqrecords(mss, full_fmt, prefix),
        fasta, "fasta")
    return fasta


def get_adam_ms_variations(unwrapper, padding):
    try:
        return AdamMSVariations.objects.get(
            unwrapper=unwrapper,
            padding=padding,
        )
    except AdamMSVariations.DoesNotExist:
        index_dir = get_unique_path()
        fasta = _build_ms_variations(unwrapper, padding)
        os.mkdir(index_dir)
        mock_msv = AdamMSVariations(
            unwrapper=unwrapper,
            padding=padding,
            index_dump_dir=index_dir,
        )
        bowtie2build(fasta, mock_msv.index_files_prefix)
        os.unlink(fasta)
        msv, c = AdamMSVariations.objects.get_or_create(
            unwrapper=unwrapper,
            padding=padding,
            defaults={"index_dump_dir": index_dir}
        )
        if not c:
            delete_files(mock_msv)
        return msv
