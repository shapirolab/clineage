
import os
from Bio import SeqIO
import pysam
import itertools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from django.db import models
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver

from misc.utils import get_unique_path

from sequencing.runs.models import Demultiplexing
from targeted_enrichment.planning.models import UGS
from targeted_enrichment.unwrapping.models import Unwrapper
from lib_prep.workflows.models import Library, BarcodedContent
from lib_prep.multiplexes.models import Panel

from sequencing.analysis.merge import pear_with_defaults
from sequencing.analysis.index_reads import bowtie2build, bowtie2_with_defaults, create_padded_fasta


def delete_files(sender, instance, **kwargs):
    # Pass false so FileField doesn't save the model.
    for filename in getattr(instance, "files", []):
        try:
            os.unlink(filename)
        except OSError:
            raise
    for dirname in getattr(instance, "dirs", []):
        try:
            os.rmdir(dirname)
        except OSError:
            raise


class SampleReads(models.Model):
    demux = models.ForeignKey(Demultiplexing)
    barcoded_content = models.ForeignKey(BarcodedContent)
    library = models.ForeignKey(Library)
    fastq1 = models.FilePathField()
    fastq2 = models.FilePathField()

    @property
    def files(self):
        yield self.fastq1
        yield self.fastq2

    def run_merge(self, merging_scheme):
        prefix = get_unique_path()
        # TODO: check by merging scheme.
        pear_with_defaults("-f", self.fastq1,
                           "-r", self.fastq2,
                           "-o", prefix)
        mr = MergedReads.objects.create(
            demux_reads=self,
            merging_scheme=merging_scheme,
            assembled_fastq="{}.assembled.fastq".format(prefix),
            discarded_fastq="{}.discarded.fastq".format(prefix),
            unassembled_forward_fastq="{}.unassembled.forward.fastq".format(prefix),
            unassembled_reverse_fastq="{}.unassembled.reverse.fastq".format(prefix),
        )
        return mr

post_delete.connect(delete_files, SampleReads)


class MergingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()


class MergedReads(models.Model):
    demux_reads = models.ForeignKey(SampleReads)
    merging_scheme = models.ForeignKey(MergingScheme)
    # TODO: Add default path?
    # TODO: Custom FASTQ field that caches #sequences
    assembled_fastq = models.FilePathField()
    discarded_fastq = models.FilePathField()
    unassembled_forward_fastq = models.FilePathField()
    unassembled_reverse_fastq = models.FilePathField()

    @property
    def files(self):
        yield self.assembled_fastq
        yield self.discarded_fastq
        yield self.unassembled_forward_fastq
        yield self.unassembled_reverse_fastq

    def included_reads_generator(self, included_reads):
        if included_reads == 'M':
            return SeqIO.parse(self.assembled_fastq, "fastq")
        elif included_reads == 'F':
            return itertools.chain(
                SeqIO.parse(self.assembled_fastq, "fastq"),
                SeqIO.parse(self.unassembled_forward_fastq, "fastq")
            )
        else:
            raise ValueError("included_reads should be one of {}".format(ReadsIndex.INCLUDED_READS_OPTIONS))

    def create_reads_index(self, included_reads, padding):
        """
        Create an index from some of these reads. included_reads should be one
        of ReadsIndex.INCLUDED_READS_OPTIONS, and chooses which of the reads we take.
        padding controls how much to pad on each side of the reads.
        """
        index_dir = get_unique_path()
        reads = self.included_reads_generator(included_reads)
        fasta = create_padded_fasta(reads, padding)
        # TODO: clean this double code.
        os.mkdir(index_dir)
        bowtie2build(fasta, os.path.join(index_dir, ReadsIndex.INDEX_PREFIX))
        ri = ReadsIndex.objects.create(
            merged_reads=self,
            included_reads=included_reads,
            index_dump_dir=index_dir,
            padding=padding,
        )
        return ri

post_delete.connect(delete_files, MergedReads)


class ReadsIndex(models.Model):
    """
    A collection of reads used for generating a bowtie2 index for primers to be searched against
    """
    INCLUDED_READS_OPTIONS = (('M', 'Only merged'), ('F', 'Merged and unassembled_forward'),)
    INDEX_PREFIX = "index"
    merged_reads = models.ForeignKey(MergedReads)
    included_reads = models.CharField(max_length=1, choices=INCLUDED_READS_OPTIONS)
    index_dump_dir = models.FilePathField(allow_files=False, allow_folders=True)
    padding = models.IntegerField(default=5)

    def included_reads_generator(self):
        return self.merged_reads.included_reads_generator(self.included_reads)

    @property
    def index_files_prefix(self):
        return os.path.join(self.index_dump_dir, self.INDEX_PREFIX)

    @property
    def files(self):
        for i in xrange(1, 5):
            yield "{}.{}.bt2".format(self.index_files_prefix, i)
        for i in xrange(1, 3):
            yield "{}.rev.{}.bt2".format(self.index_files_prefix, i)

    @property
    def dirs(self):
        yield self.index_dump_dir

post_delete.connect(delete_files, ReadsIndex)


class UGSAssignment(models.Model):
    reads_index = models.ForeignKey(ReadsIndex)
    panel_fasta = models.FilePathField(null=True)
    primer_reads_alignment = models.FilePathField(null=True)

    @property
    def primer_reads_alignment_name(self):
        return 'primer_reads_alignment_{}.sam'.format(self.id)

    @property
    def unwrappers(self):
        return self.reads_index.merged_reads.demux_reads.library.unwrappers

    @staticmethod
    def left_unwrapper_name(unwrapper):
        return 'unwrapper_{}_left'.format(unwrapper.id)

    @staticmethod
    def right_unwrapper_name(unwrapper):
        return 'unwrapper_{}_right'.format(unwrapper.id)

    @staticmethod
    def trim_one_base_from_left(seq_recs):
        """
        Workaround: It turns out that Adam was generating the primer sequences using "bedtools getfasta" which has a
        different interpretation of indices to that of UCSD and our implementation.
        """
        for sr in seq_recs:
            yield sr[1:]

    def primers_seqrecords_generator(self):
        for uw in self.unwrappers:
            yield SeqRecord(Seq(uw.left_margin.seq), id=self.left_unwrapper_name(uw), name='', description='')
            yield SeqRecord(Seq(uw.right_margin.seq), id=self.right_unwrapper_name(uw), name='', description='')
            # name='', description='' are workarounds for the '<unknown description>' that is being outputted and
            # breaks the downstream bowtie2 alignment.

    def create_panel_fasta(self):
        panel_fasta_name = get_unique_path("fasta")
        SeqIO.write(self.trim_one_base_from_left(self.primers_seqrecords_generator()), panel_fasta_name, "fasta")
        self.panel_fasta = panel_fasta_name
        self.save()

    def align_primers_to_reads(self):
        primer_reads_alignment_name = get_unique_path("sam")
        bowtie2_with_defaults('-x', self.reads_index.index_files_prefix,
                              '-f', self.panel_fasta,
                              '-S', primer_reads_alignment_name)
        self.primer_reads_alignment = primer_reads_alignment_name
        self.save()

    def read_sam(self):
        with pysam.AlignmentFile(self.primer_reads_alignment, "rb") as samfile:
            for r in samfile:
                if r.is_unmapped:
                    continue  # unmapped TODO: dump in appropriate bin
                read_name = samfile.getrname(r.reference_id)
                primer_name = r.query_name
                yield read_name, primer_name

    def collect_mappings_from_sam(self):
        reads_by_primers = {}
        for read_name, primer_name in self.read_sam():
            reads_by_primers.setdefault(primer_name, set()).add(read_name)
        return reads_by_primers

    def read_ids_by_unwrapper(self):
        reads_by_ugs = self.collect_mappings_from_sam()
        for uw in self.unwrappers:
            # TODO: shirgun
            yield uw, set(reads_by_ugs.get(self.left_unwrapper_name(uw), [])) & set(reads_by_ugs.get(self.right_unwrapper_name(uw), []))

    def reads_by_unwrapper(self):
        fq_dict = SeqIO.to_dict(self.reads_index.included_reads_generator())

        for uw, read_ids in self.read_ids_by_unwrapper():
            reads = [fq_dict[read_id] for read_id in read_ids]
            yield uw, reads


class SequencingData(models.Model):  # This contains the actual data.
    merged_reads = models.ForeignKey(MergedReads)
    unwrapper = models.ForeignKey(Unwrapper)
    target_offset = models.IntegerField(null=True)
    fastq = models.FilePathField(null=True)
    vcf = models.FilePathField(null=True)

    class Meta:
        index_together = (
            ("merged_reads", "unwrapper"),
        )
