
from plumbum import local
import csv
import io
import os
import resource
from Bio import SeqIO

from misc.utils import get_unique_path, unique_dir_cm, unlink

from django.db.models.query import QuerySet

from sequencing.runs.models import Demultiplexing
from sequencing.analysis.models import SampleReads


def run_bcl2fastq(bcl_folder, sample_sheet_path, fastq_folder):
    assert resource.getrlimit(resource.RLIMIT_OFILE)[1] == 4096
    resource.setrlimit(resource.RLIMIT_OFILE, (4096, 4096))
    bcl2fastq = local["bcl2fastq"]
    bcl2fastq_with_defaults = bcl2fastq["--no-lane-splitting",
                                        "--processing-threads", 20]
    bcl2fastq_with_defaults(
        "--runfolder-dir", bcl_folder,
        "--sample-sheet", sample_sheet_path,
        "--output-dir", fastq_folder,
    )


SAMPLE_ID = "Sample_ID"
SAMPLE_NAME = "Sample_Name"
RIGHT_BARCODE_ID = "I7_Index_ID"
RIGHT_BARCODE_SEQ = "index"
LEFT_BARCODE_ID = "I5_Index_ID"
LEFT_BARCODE_SEQ = "index2"

SAMPLESHEET_HEADERS = [
    SAMPLE_ID,
    SAMPLE_NAME,
    "Sample_Plate",
    "Sample_Well",
    RIGHT_BARCODE_ID,
    RIGHT_BARCODE_SEQ,
    LEFT_BARCODE_ID,
    LEFT_BARCODE_SEQ,
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

def generate_sample_sheets(bcs, name, date, description, kit, demux_scheme):
    read_length = kit.read_length
    fwd_read_adaptor = kit.fwd_read_adaptor
    rev_read_adaptor = kit.rev_read_adaptor
    sio = io.StringIO()
    sio.write(HEADER_FORMAT.format(
        run_name=name,
        run_date=date,
        run_desc=description,
        read_length=read_length,
        fwd_read_adaptor=fwd_read_adaptor,
        rev_read_adaptor=rev_read_adaptor,
    ))
    dia = csv.unix_dialect()
    dia.quoting = csv.QUOTE_NONE
    w = csv.DictWriter(sio, fieldnames=SAMPLESHEET_HEADERS, dialect=dia)
    w.writeheader()
    for bc in bcs:
        w.writerow({
            SAMPLE_ID: bc.id,
            SAMPLE_NAME: bc.id,
            LEFT_BARCODE_ID: bc.barcodes.left.name,
            LEFT_BARCODE_SEQ: bc.barcodes.left.ref_sequence,
            RIGHT_BARCODE_ID: bc.barcodes.right.name,
            RIGHT_BARCODE_SEQ: bc.barcodes.right.ref_sequence,
        })
    s = sio.getvalue()
    sio.close()
    return s


UNDETERMINED = "Undetermined"
SAMPLE_FASTQ_GZ_FORMAT = "{id}_S{idx}_R{read:d}_001.fastq.gz"
SAMPLE_FASTQ_FORMAT = "{id}_S{idx}_R{read:d}_001.fastq"


def run_demux(ngs_run, demux_scheme):
    bc_l = []
    bc_libs_d = {}
    bc_idx_d = {}
    for library in ngs_run.libraries.select_subclasses():
        it = library.barcoded_contents
        if isinstance(it, QuerySet):
            it = it.select_related(
                "barcodes",
                "barcodes__left",
                "barcodes__right"
            )
        for bc in it:
            bc_libs_d[bc] = library
            bc_l.append(bc)
            bc_idx_d[bc] = len(bc_l)  # NOTE: this off-by-1 is intentional.
    description = "_".join(lib.name for lib in ngs_run.libraries.all())
    sample_sheet = generate_sample_sheets(
        bcs=bc_l,
        name=ngs_run.name,
        date=ngs_run.date,
        description=description,
        kit=ngs_run.kit,
        demux_scheme=demux_scheme,
    )
    with unique_dir_cm() as fastq_folder:
        with unlink(get_unique_path("csv")) as sample_sheet_path:
            with open(sample_sheet_path, "w") as f:
                f.write(sample_sheet)
            run_bcl2fastq(ngs_run.bcl_directory, sample_sheet_path, fastq_folder)
        demux = Demultiplexing.objects.create(
            ngs_run=ngs_run,
            demux_scheme=demux_scheme,
        )
        files = []
        # TODO: nicer handling
        for read in [1, 2]:
            os.unlink(os.path.join(fastq_folder, SAMPLE_FASTQ_GZ_FORMAT.format(
                id=UNDETERMINED,
                idx=0,
                read=read,
            )))
        for bc, idx in bc_idx_d.items():
            for fastqgz in [os.path.join(fastq_folder, 
                SAMPLE_FASTQ_GZ_FORMAT.format(
                    id=bc.id,
                    idx=idx,
                    read=read,
                )) for read in [1, 2]]:
                local["gunzip"](fastqgz)
            fastq1, fastq2 = [os.path.join(fastq_folder, 
                SAMPLE_FASTQ_FORMAT.format(
                    id=bc.id,
                    idx=idx,
                    read=read,
                )) for read in [1, 2]]
            files.append(fastq1)
            files.append(fastq2)
            with open(fastq1) as f:
                seq = SeqIO.parse(f, format="fastq")
                num_reads = sum(1 for x in seq)
            yield SampleReads.objects.create(
                demux=demux,
                barcoded_content=bc,
                library=bc_libs_d[bc],
                num_reads=num_reads,
                fastq1=fastq1,
                fastq2=fastq2,
            )
