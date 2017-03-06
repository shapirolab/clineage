# __author__ = 'veronika'

import os
import uuid
import plumbum
import vcf
import itertools
import contextlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from distributed import as_completed
from Bio.Sequencing.Applications import BwaAlignCommandline, BwaSampeCommandline, \
     SamtoolsViewCommandline, SamtoolsVersion1xSortCommandline, BwaIndexCommandline, SamtoolsFaidxCommandline

from misc.utils import get_unique_path, unlink, unique_file_cm, unique_dir_cm, get_get_or_create
from django.contrib.auth.models import User

from sequencing.analysis.snps.models import ReadsAlignment, VCFReads, SNPReads
from sequencing.analysis.snps.vcf_generator import mpileup, bcf_tools
from sequencing.analysis.snps.parse_snps import convert_snp_dict_to_table, create_snp_file
from sequencing.analysis.snps.models import _parse_vcf_file, AmpliconCollectionBWAIndex
from sequencing.analysis.models_common import BWAIndexMixin, SampleReads
from targeted_enrichment.amplicons.models import AmpliconCollection, Amplicon
from targeted_enrichment.planning.models import SNP

rm = plumbum.local["rm"]


def sort_and_index_bam_file(bam_file):
    """
    Sort and index BAM file
    :param bam_file: bam file
    :return: sorted bam file
    """

    with unique_dir_cm() as samtools_sort_data:
        with unique_file_cm(ext='bam') as bam_sort_file:
            samtools_sort_cmd = SamtoolsVersion1xSortCommandline(input=bam_file,
                                                                 T=samtools_sort_data,
                                                                 o=bam_sort_file,
                                                                 O='bam')
            stdout, stderr = samtools_sort_cmd()
            unlink(samtools_sort_data)
            os.unlink(bam_file)
            return bam_sort_file


def sam_to_bam(sam_file):
    """
    Convert SAM to BAM
    :param sam_file: Sam file
    :return: Bam file
    """
    with unique_file_cm(ext='bam') as bam_file:
        samtools_view_cmd = SamtoolsViewCommandline(input_file=sam_file, S=True, b=True, u=True, o=bam_file)
        samtools_view_cmd()

        rm(sam_file)
        return bam_file


def snp_alignment(sample_read, ref_fasta):
    """
    Using samtools to align reads from paired fastq to reference seq
    :param sample_read: SampleReads object
    :param ref_fasta: path to fasta and index of reference sequence
    :return: bam file
    """

    with unique_file_cm(ext='sai') as output_sai_file1:
        with unique_file_cm(ext='sai') as output_sai_file2:
            output_sai_command1 = BwaAlignCommandline(reference=ref_fasta, read_file=sample_read.fastq1)
            stdout, stderr = output_sai_command1(stdout=output_sai_file1)
            print('Output R1 ', stderr)

            output_sai_command2 = BwaAlignCommandline(reference=ref_fasta, read_file=sample_read.fastq2)
            stdout, stderr = output_sai_command2(stdout=output_sai_file2)
            print('Output R2 ', stderr)

            sampe_cmd = BwaSampeCommandline(reference=ref_fasta, sai_file1=output_sai_file1, sai_file2=output_sai_file2,
                                            read_file1=sample_read.fastq1, read_file2=sample_read.fastq2,
                                            )
            with unique_file_cm(ext='sam') as sampe_data:
                sampe_cmd(stdout=sampe_data)
                bam_file = sam_to_bam(sampe_data)

                #  TODO: remove
                assert not os.path.exists(sampe_data)
                assert os.path.exists(output_sai_file1)
                assert os.path.exists(output_sai_file2)
                rm(output_sai_file1)
                rm(output_sai_file2)
                assert not os.path.exists(output_sai_file1)
                assert not os.path.exists(output_sai_file2)

                return bam_file


def get_amplicon_seqrecord(seq_fmt, prefix):
        # name='', description='' are workarounds for the '<unknown
        # description>' that is being outputted otherwise
        return SeqRecord(
            Seq(seq_fmt),
            id="{}".format(prefix),
            name='',
            description=''
        )


def make_faidx(fasta_file):
    """
    makes faidx in the fasta folder
    :param fasta_file: fasta file for indexing
    :return: faidx file
    """

    samtools_faidx_cmd = SamtoolsFaidxCommandline(reference=fasta_file)
    samtools_faidx_cmd()
    faidx_file = '{}.fai'.format(fasta_file)
    return faidx_file


def prepare_seq_for_fasta(amplicon_collection):
    """
    Template for fasta lines
    :param amplicon_collection: amplicons collection object
    :return: lines for fasta file
    """
    for amp in amplicon_collection.amplicons.all():
        fmt = amp.slice.sequence
        if amp.subclass.left_ugs.slice.start_pos == amp.slice.start_pos:
            full_fmt = "NNN{slice}NNN".format(
                slice=fmt,
            )
        else:
            full_fmt = "{left}{slice}{right}".format(
                left=amp.subclass.left_margin,
                slice=fmt,
                right=amp.subclass.right_margin,
            )
        prefix = "{}".format(amp.id)
        yield get_amplicon_seqrecord(full_fmt, prefix)


@contextlib.contextmanager
def get_fasta_file(amplicon_collection, path=None):
    """
    Creates fasta file out of amplicons collection
    :param amplicon_collection:
    :param path:
    :return: fasta file
    """
    try:
        fasta = str(uuid.uuid4())
        fasta += ".{}".format("fa")
        if path is not None:
            fasta = os.path.join(path, fasta)
        else:
            fasta = get_unique_path("fa")
        SeqIO.write(prepare_seq_for_fasta(amplicon_collection), fasta, "fasta")
        yield fasta
    except:
        try:
            os.unlink(fasta)
        except:
            pass
        raise


def fiter_amplicons_for_snp_analysis(amplicons):
    """
    Filter amplicons for SNP amplicons only
    :param amplicon: amplicon list
    :return: list of filtered amplicons
    """
    amp_for_fasta = []
    for amp in amplicons:
        for s in amp.slice.contains.all():
            for t in s.target_set.all():
                try:
                    snp = t.snp
                except SNP.DoesNotExist:
                    continue
                amp_for_fasta.append(snp)
    return amp_for_fasta


def get_snp_amplicon_collection(amplicons):
    """
    Create amplicons collection object for the snp run
    :param amplicons: amplicons list
    :return: AmpliconCollection object
    """
    amplicon_collection = amplicons.amplicons.all()

    def inner(raise_or_create_with_defaults):
        amplicon_collection = fiter_amplicons_for_snp_analysis(amplicons)

        return raise_or_create_with_defaults(
            amplicons=amplicon_collection,
        )

    return get_get_or_create(inner, AmpliconCollection,
                             amplicons=amplicon_collection,
                             )


def check_targets_id(vcf_file):
    cell_snps = []
    snps_rec = []
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for r in vcf_reader:
        if len(r.REF) == 1 and len(r.ALT) != 1:
            snps_rec.append(r)
    for snp in snps_rec:
        amp = Amplicon.objects.get(id=snp.CHROM)
        cell_snps.append(amp.id)

    return cell_snps


def get_alignment_reference(amplicon_collection):
    """
    Creates index for amplicon_collection and creates AmpliconCollectionBWAIndex object for the index
    :param amplicon_collection: amplicons collection object
    :return: AmpliconCollectionBWAIndex object
    """
    def inner(raise_or_create_with_defaults):
        with unique_dir_cm() as index_dir:
            bwa_index = BWAIndexMixin(index_dump_dir=index_dir)  # hold IndexMixin later used in defaults

            with get_fasta_file(amplicon_collection, path=index_dir) as fasta:
                index_cmd = BwaIndexCommandline(infile=fasta, prefix=bwa_index.index_files_prefix, algorithm="bwtsw")
                index_cmd()
                faidx_file = make_faidx(fasta)
                assert os.path.dirname(fasta) == index_dir
                bwa_index.index_dump_dir = os.path.dirname(fasta)
                assert os.path.commonprefix([fasta, bwa_index.index_dump_dir]) == bwa_index.index_dump_dir
                assert os.path.exists(fasta)
                return raise_or_create_with_defaults(
                    index_dump_dir=index_dir,
                    fasta_file=fasta,
                    faidx_file=faidx_file,
                )

    return get_get_or_create(inner, AmpliconCollectionBWAIndex,
                             amplicon_collection=amplicon_collection,
                             )


def align_reads(sample_read, alignment_reference):
    """
    Aligning reads with ref fasta, sorting the BAM file and creating ReadsAlignment object
    :param sample_read: SampleReads object
    :param alignment_reference: reference fasta
    :return: ReadsAlignment object
    """

    def inner(raise_or_create_with_defaults):
        alignment_file = snp_alignment(sample_read, ref_fasta=alignment_reference.index_files_prefix)  # TODO: get_get_or_create
        sort_bam_file = sort_and_index_bam_file(alignment_file)

        # os.unlink(alignment_file)  # TODO: try unlink
        assert not os.path.exists(alignment_file)
        return raise_or_create_with_defaults(
            bam_file=sort_bam_file,
        )

    return get_get_or_create(inner, ReadsAlignment,
                             sample_read=sample_read,
                             alignment_reference=alignment_reference,
                             )


def pileup(reads_alignment):
    """
    Applying mpileup function on BAM files and creating vcf object.
    :param reads_alignment: ReadsAlignment object
    :return: vcf_model object
    """

    def inner(raise_or_create_with_defaults):
        # converting SNP reads to VCF  files
        vcf_file = mpileup(reads_alignment.bam_file, reads_alignment.alignment_reference.fasta_file)

        return raise_or_create_with_defaults(
            vcf_file=vcf_file,
        )

    return get_get_or_create(inner, VCFReads,
                             reads_alignment=reads_alignment,
                             )


def snp_table(cell_snps_list, score_threshold=0.5, method='explicit_snps', min_cover=5, partner=None):
    snp_dict = {}
    if score_threshold < 0.5 or score_threshold > 1:
        score_threshold = 0.5
        print('WARN: scorethreshold is floored to 0.5')

    if partner:
        partnet_id = User.objects.get(username=partner)

    for cell_snps in cell_snps_list:

        cell_snps_dict = cell_snps.extract_snp_file()
        cell_id = cell_snps.vcf_read.reads_alignment.sample_read.id
        if partner:
            sr = SampleReads.objects.get(id=cell_id)
            if sr.barcoded_content.subclass.content.cell.individual.partner != partnet_id:
                continue

        for chrom, loc in cell_snps_dict:
            if (chrom, loc) not in snp_dict:
                snp_dict[(chrom, loc)] = {}
            snp_dict[(chrom, loc)].setdefault(str(cell_id), cell_snps_dict[(chrom, loc)])

    rows = convert_snp_dict_to_table(snp_dict,
                                     score_threshold=score_threshold,
                                     method=method,
                                     min_cover=min_cover)
    return snp_dict, rows


def parse_vcf(vcf_obj, min_cover=1):
    """
    Creating SNP's object out of VCF file
    :param vcf_obj: VCFReads object
    :param min_cover: THe minimum cover of reads per SNP
    :return: SNPReads obj
    """

    def inner(raise_or_create_with_defaults):
        # converting VCF file to SNP dict and saving as pickle
        snps_dict = _parse_vcf_file(vcf_obj.vcf_file, min_cover=min_cover)
        snps_dict_file = create_snp_file(snps_dict)

        return raise_or_create_with_defaults(
            snps_dict=snps_dict_file,
        )

    # return (cell_snps, vcf_obj.reads_alignment.sample_read.id)
    return get_get_or_create(inner, SNPReads,
                             vcf_read=vcf_obj,
                             min_cover=min_cover,
                             )

################# Parallel ##############################################


def run_parallel(executor, sample_reads, amplicon_collection, min_cover=5):
    # *currently in dworker.q

    alignment_reference = executor.submit(get_alignment_reference, amplicon_collection)
    read_alignment = executor.map(align_reads, sample_reads, itertools.repeat(alignment_reference))
    yield read_alignment

    vcf_list = executor.map(pileup, read_alignment)
    yield vcf_list

    snp_list = executor.map(parse_vcf, vcf_list, min_cover=min_cover)
    yield snp_list
