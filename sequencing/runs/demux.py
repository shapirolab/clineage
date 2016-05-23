
from plumbum import local
import csv
import io
import os

from misc.utils import get_unique_path


def run_bcl2fastq(bcl_folder, sample_sheet_path):
    bcl2fastq = local["bcl2fastq"]
    bcl2fastq_with_defaults = bcl2fastq["--no-lane-splitting",
                                            "--processing-threads", 20]
    bcl2fastq_with_defaults(
        "--runfolder-dir", bcl_folder,
        "--sample-sheet", sample_sheet_path)


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
    w = csv.DictWriter(sio, fieldnames=SAMPLESHEET_HEADERS)
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
SAMPLE_FASTQ_GZ_FORMAT = "{id}_S{idx:03d}_L001_R{read:d}_001.fastq.gz"
SAMPLE_FASTQ_FORMAT = "{id}_S{idx:03d}_L001_R{read:d}_001.fastq"


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
    description = "_".join(lib.name for lib in ngsrun.libraries.all())
    sample_sheet = generate_sample_sheets(
        bcs=bc_l,
        name=ngsrun.name,
        date=ngsrun.date,
        description=description,
        kit=ngsrun.kit,
        demux_scheme=demux_scheme,
        max_samples=max_samples
    )
    sample_sheet_path = get_unique_path("csv")
    with open(sample_sheet_path, "w") as f:
        f.write(sample_sheet)
    fastq_folder = get_unique_path()
    os.mkdir(fastq_folder)
    run_bcl2fastq(fastq_folder, sample_sheet_path)
    demux = Demultiplexing.objects.create(
        ngs_run=ngs_run,
        demux_scheme=demux_scheme,
    )
    files = []
    # TODO: nicer handling
    for read in [1, 2]:
        os.unlink(os.path.join(sample_sheet_path, SAMPLE_FASTQ_GZ_FORMAT.format(
            id=UNDETERMINED,
            idx=0,
            read=read,
        )))
    for bc, idx in bc_idx_d.items():
        for fastqgz in [os.path.join(sample_sheet_path, 
            SAMPLE_FASTQ_GZ_FORMAT.format(
                id=bc.id,
                idx=idx,
                read=read,
            )) for read in [1, 2]]:
            local["gunzip"](fastqgz)
        fastq1, fastq2 = [os.path.join(sample_sheet_path, 
            SAMPLE_FASTQ_FORMAT.format(
                id=bc.id,
                idx=idx,
                read=read,
            )) for read in [1, 2]]
        files.append(fastq1)
        files.append(fastq2)
        yield SampleReads.objects.create(
            demux=demux,
            barcoded_content=bc,
            library=bc_libs_d[bc],
            fastq1=fastq1,
            fastq2=fastq2,
        )
