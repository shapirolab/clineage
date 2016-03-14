import os

from sequencing.analysis.models import SampleReads, AdamMergedReads, \
    AdamReadsIndex, AdamMarginAssignment, AdamAmpliconReads, \
    AdamMSVariations, AdamHistogram
from misc.utils import get_unique_path

file_fixtures_path = os.path.join(*(os.path.split(os.path.dirname(os.path.realpath(__file__)))[:-1] + ("ngs_fixtures",)))
adam_fixtures_path = os.path.join(file_fixtures_path, "adam")
adam_bcs_fixtures_path = os.path.join(adam_fixtures_path, "bcs")
adam_mss_fixtures_path = os.path.join(adam_fixtures_path, "mss")


def generate_samplereads(request, bc):
    fastq_r1 = get_unique_path("fastq")
    fastq_r2 = get_unique_path("fastq")
    os.symlink(
        os.path.join(adam_bcs_fixtures_path, bc, "{}_R1.fastq".format(bc)),
        fastq_r1
    )
    os.symlink(
        os.path.join(adam_bcs_fixtures_path, bc, "{}_R2.fastq".format(bc)),
        fastq_r2
    )
    sr = SampleReads.objects.create(
        demux=request.getfuncargvalue("demultiplexing"),
        barcoded_content=request.getfuncargvalue("magicalpcr1barcodedcontent"),
        library=request.getfuncargvalue("magicalpcr1library"),
        fastq1=fastq_r1,
        fastq2=fastq_r2,
    )
    return sr


def generate_mergedreads(request, bc):
    src_prefix = os.path.join(adam_bcs_fixtures_path, bc, bc)
    dst_prefix = get_unique_path()
    os.symlink(
        "{}.assembled.fastq".format(src_prefix),
        "{}.assembled.fastq".format(dst_prefix)
    )
    os.symlink(
        "{}.discarded.fastq".format(src_prefix),
        "{}.discarded.fastq".format(dst_prefix)
    )
    os.symlink(
        "{}.unassembled.forward.fastq".format(src_prefix),
        "{}.unassembled.forward.fastq".format(dst_prefix)
    )
    os.symlink(
        "{}.unassembled.reverse.fastq".format(src_prefix),
        "{}.unassembled.reverse.fastq".format(dst_prefix)
    )
    mr = AdamMergedReads.objects.create(
        sample_reads=generate_samplereads(request, bc),
        assembled_fastq="{}.assembled.fastq".format(dst_prefix),
        discarded_fastq="{}.discarded.fastq".format(dst_prefix),
        unassembled_forward_fastq="{}.unassembled.forward.fastq".format(dst_prefix),
        unassembled_reverse_fastq="{}.unassembled.reverse.fastq".format(dst_prefix),
    )
    return mr


def link_index_files(src_dir, dst_dir):
    os.mkdir(dst_dir)
    for index_file in ["index.{}.bt2".format(s) for s in \
            ["1","2","3","4","rev.1","rev.2"]]:
        os.symlink(
            os.path.join(src_dir, index_file),
            os.path.join(dst_dir, index_file),
        )


def generate_readsindex(request, bc, included_reads):
    src_dir = os.path.join(adam_bcs_fixtures_path, bc, included_reads)
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    ri = AdamReadsIndex.objects.create(
        merged_reads=generate_mergedreads(request, bc),
        included_reads=included_reads,  # Only merged
        index_dump_dir=dst_dir,
        padding=5,
    )
    return ri


def generate_adammarginassignment(request, bc, inc):
    request.getfuncargvalue("require_amplicons")
    alignment_file_name = get_unique_path("sam")
    os.symlink(
        os.path.join(adam_bcs_fixtures_path, bc, inc, "margins.sam"),
        alignment_file_name
    )
    ama = AdamMarginAssignment.objects.create(
        reads_index=generate_readsindex(request, bc, inc),
        assignment_sam=alignment_file_name,
    )
    return ama


def generate_adamampliconreads(request, bc, inc, unw):
    fastq_path = get_unique_path("fastq")
    base_fixtures_path = os.path.join(adam_bcs_fixtures_path, bc, inc, unw, unw)
    os.symlink(
        "{}.fastq".format(base_fixtures_path),
        fastq_path
    )
    fastq1_path = get_unique_path("fastq")
    os.symlink(
        "{}_R1.fastq".format(base_fixtures_path),
        fastq1_path
    )
    fastq2_path = get_unique_path("fastq")
    os.symlink(
        "{}_R2.fastq".format(base_fixtures_path),
        fastq2_path
    )
    aar = AdamAmpliconReads.objects.create(
        margin_assignment=generate_adammarginassignment(request, bc, inc),
        amplicon=request.getfuncargvalue("pu_{}".format(unw)),
        fastq=fastq_path,
        fastq1=fastq1_path,
        fastq2=fastq2_path,
    )
    return aar


def generate_amsv(request, unw):
    src_dir = os.path.join(adam_mss_fixtures_path, unw)
    dst_dir = get_unique_path()
    link_index_files(src_dir, dst_dir)
    amsv = AdamMSVariations.objects.create(
        amplicon=request.getfuncargvalue("pu_{}".format(unw)),
        index_dump_dir=dst_dir,
        padding=50,
    )
    return amsv


def generate_adamhistogram(request, bc, inc, unw):
    alignment_file_name = get_unique_path("sam")
    os.symlink(
        os.path.join(adam_bcs_fixtures_path, bc, inc, unw,
            "{}.sam".format(unw)),
        alignment_file_name
    )
    adamampliconreads = generate_adamampliconreads(request, bc, inc, unw)
    ama = AdamHistogram.objects.create(
        sample_reads=adamampliconreads.margin_assignment.reads_index \
            .merged_reads.sample_reads,
        amplicon_reads=adamampliconreads,
        assignment_sam=alignment_file_name,
    )
    return ama
