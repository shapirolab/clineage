
import itertools
import csv
import io

from django.contrib.auth.models import User
from django.db import models
from django.db.models.query import QuerySet
from lib_prep.workflows.models import Library, BarcodedContent
from primers.parts.models import IlluminaReadingAdaptor1, \
    IlluminaReadingAdaptor2

SAMPLESHEET_HEADERS = [
"Sample_ID",
"Sample_Name",
"Sample_Plate",
"Sample_Well",
"I7_Index_ID",
"index",
"I5_Index_ID",
"index2",
"Sample_Project",
"Description",
]

HEADER_FORMAT = \
"""[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Experiment Name,{run_name},,,,,,,,
Date,{run_date},,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,{run_desc},,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
{read_length},,,,,,,,,
{read_length},,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,{fwd_read_adaptor},,,,,,,,
AdapterRead2,{rev_read_adaptor},,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
"""
class MachineType(models.Model):
    company = models.CharField(max_length=50)
    model = models.CharField(max_length=50)

    def __unicode__(self):
        return self.company + '_' + self.model


class Machine(models.Model):
    machineid = models.CharField(max_length=50)
    name = models.CharField(max_length=100, null=True, blank=True)
    type = models.ForeignKey(MachineType)

    def __unicode__(self):
        return self.type.__unicode__() + '_' + self.machineid


class NGSKit(models.Model):
    reading_adaptor1 = models.ForeignKey(IlluminaReadingAdaptor1)
    reading_adaptor2 = models.ForeignKey(IlluminaReadingAdaptor2)
    read_length = models.IntegerField(null=True)

    @property
    def fwd_read_adaptor(self):
        """
        Return the value for the Adapter field in the SampleSheet.
        """
        return self.reading_adaptor2.sequence[1:].rev_comp()

    @property
    def rev_read_adaptor(self):
        """
        Return the value for the AdapterRead2 field in the SampleSheet.
        """
        return self.reading_adaptor1.sequence.rev_comp()


class NGSRun(models.Model):
    libraries = models.ManyToManyField(Library)
    bcl_directory = models.FilePathField(allow_files=False, allow_folders=True, null=True)
    name = models.CharField(max_length=100, unique=True)
    machine = models.ForeignKey(Machine)
    kit = models.ForeignKey(NGSKit, null=True)
    user = models.ForeignKey(User)
    date = models.DateField()

    def generate_sample_sheets(self, demux_scheme, max_samples=None):
        """
        Create samplesheet(s) for this given run for demultiplexing with given
        demux_scheme. If max_samples is specified, creates multiple
        samplesheets, each with at most max_samples samples. Otherwise creates
        a single samplesheet.
        Currently, demux_scheme doesn't change anything.
        FIXME: Use demux_scheme.
        """
        out = []

        def rows_iter():
            for library in self.libraries.select_subclasses():
                it = library.barcoded_contents
                if isinstance(it, QuerySet):
                    it = it.select_related(
                        "barcodes",
                        "barcodes__left",
                        "barcodes__right"
                    )
                for bc in library.barcoded_contents:
                    yield {
                        "Sample_ID": bc.id,
                        "Sample_Name": bc.id,
                        "I7_Index_ID": bc.barcodes.left.name,
                        "index": bc.barcodes.left.sequence,
                        "I5_Index_ID": bc.barcodes.right.name,
                        "index2": bc.barcodes.right.sequence,  # Or ref_seq?
                    }
        rows = rows_iter()
        while True:
            bio = io.BytesIO()
            bio.write(HEADER_FORMAT.format(
                run_name=self.name,
                run_date=self.date,
                run_desc="_".join(lib.name for lib in self.libraries.all()),
                read_length=self.kit.read_length,
                fwd_read_adaptor=self.kit.fwd_read_adaptor,
                rev_read_adaptor=self.kit.rev_read_adaptor,
            ))
            w = csv.DictWriter(bio, fieldnames=SAMPLESHEET_HEADERS)
            w.writeheader()
            if max_samples is None:
                for row in rows:
                    w.writerow(row)
                b = bio.getvalue()
                bio.close()
                return b
            else:
                row = None
                for row in itertools.islice(rows, max_samples):
                    w.writerow(row)
                if row is None:
                    return out
            out.append(bio.getvalue())
            bio.close()
