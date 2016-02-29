import os
from Bio import SeqIO
import pysam
from itertools import chain
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import uuid

from django.db import models
from django.conf import settings

from sequencing.runs.models import Demultiplexing
from targeted_enrichment.planning.models import UGS
from targeted_enrichment.unwrapping.models import Unwrapper
from lib_prep.workflows.models import Library, BarcodedContent
from lib_prep.multiplexes.models import Panel

from sequencing.analysis.merge import pear_with_defaults
from sequencing.analysis.index_reads import bowtie2build, bowtie2_with_defaults


class DemultiplexedReads(models.Model):
    demux = models.ForeignKey(Demultiplexing)
    barcoded_content = models.ForeignKey(BarcodedContent)
    library = models.ForeignKey(Library)
    fastq1 = models.FilePathField()
    fastq2 = models.FilePathField()

    def run_merge(self, merging_scheme):
        prefix = os.path.join(settings.DATA_STORE,"{}".format(uuid.uuid4()))
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


class MergingScheme(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()


class MergedReads(models.Model):
    demux_reads = models.ForeignKey(DemultiplexedReads)
    merging_scheme = models.ForeignKey(MergingScheme)
    # TODO: Add default path?
    # TODO: Custom FASTQ field that caches #sequences
    assembled_fastq = models.FilePathField()
    discarded_fastq = models.FilePathField()
    unassembled_forward_fastq = models.FilePathField()
    unassembled_reverse_fastq = models.FilePathField()


class ReadsIndex(models.Model):
    """
    A collection of reads used for generating a bowtie2 index for primers to be searched against
    """
    INCLUDED_READS_OPTIONS = (('M', 'Only merged'), ('F', 'Merged and unassembled_forward'),)
    merged_reads = models.ForeignKey(MergedReads)
    included_reads = models.CharField(max_length=1, choices=INCLUDED_READS_OPTIONS)
    index_dump_dir = models.FilePathField(allow_files=False, allow_folders=True, null=True)
    padding = models.IntegerField(default=5)

    @property
    def padded_reads_fasta(self):
        return os.path.join(self.index_dump_dir,"reads.fasta")
    
    @property
    def index_files_prefix(self):
        return os.path.join(self.index_dump_dir,"index")

    def pad_records(self, records):
        for record in records:
            yield "N"*self.padding + record + "N"*self.padding

    def included_reads_generator(self):
        if self.included_reads == 'M':
            return SeqIO.parse(self.merged_reads.assembled_fastq, "fastq")
        if self.included_reads == 'F':
            return chain(
                SeqIO.parse(self.merged_reads.assembled_fastq, "fastq"),
                SeqIO.parse(self.merged_reads.unassembled_forward_fastq, "fastq")
            )
        #IntegrityError
        raise Exception()

    def create_final_merged_fastq(self):
        file_path = os.path.join(settings.DATA_STORE, "{}".format(uuid.uuid4()))
        os.mkdir(file_path)
        padded_reads = self.pad_records(self.included_reads_generator())
        self.index_dump_dir = file_path
        SeqIO.write(padded_reads, self.padded_reads_fasta, "fasta")
        self.save()

    def collect_bt_files(self):
        for path in os.listdir(self.index_dump_dir):
            if path[:len(self.index_files_prefix)] != self.index_files_prefix:
                continue
            yield path

    def index_reads(self):
        bowtie2build(self.padded_reads_fasta, self.index_files_prefix)


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
        panel_fasta_name = os.path.join(settings.DATA_STORE, "{}.fasta".format(format(uuid.uuid4())))
        SeqIO.write(self.trim_one_base_from_left(self.primers_seqrecords_generator()), panel_fasta_name, "fasta")
        self.panel_fasta = panel_fasta_name
        self.save()

    def align_primers_to_reads(self):
        primer_reads_alignment_name = os.path.join(settings.DATA_STORE, "{}.sam".format(format(uuid.uuid4())))
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


class SequencingData(models.Model): # This contains the actual data.
    merged_reads = models.ForeignKey(MergedReads)
    unwrapper = models.ForeignKey(Unwrapper)
    target_offset = models.IntegerField(null=True)
    fastq = models.FilePathField(null=True)
    vcf = models.FilePathField(null=True)

    class Meta:
        index_together = (
            ("merged_reads", "unwrapper"),
        )
