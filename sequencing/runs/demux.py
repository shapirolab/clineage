
from plumbum import local
import itertools
import csv
import io


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


def run_bcl2fastq(bcl_folder, sample_sheet_path):
    bcl2fastq = local["bcl2fastq"]
    bcl2fastq_with_defaults = bcl2fastq["--no-lane-splitting",
                                            "--processing-threads", 20]
    bcl2fastq_with_defaults(
        "--runfolder-dir", bcl_folder,
        "--sample-sheet", sample_sheet_path)

def _rows_iter(bcs):
    for bc in bcs:
        yield {
            "Sample_ID": bc.id,
            "Sample_Name": bc.id,
            "I7_Index_ID": bc.barcodes.left.name,
            "index": bc.barcodes.left.sequence,
            "I5_Index_ID": bc.barcodes.right.name,
            "index2": bc.barcodes.right.sequence,  # Or ref_seq?
        }

def generate_sample_sheets(barcodes, name, date, description, read_length, fwd_read_adaptor, rev_read_adaptor, demux_scheme, max_samples=None):
        out = []
        rows = _rows_iter(barcodes)
        while True:
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
            if max_samples is None:
                for row in rows:
                    w.writerow(row)
                b = sio.getvalue()
                sio.close()
                return b
            else:
                row = None
                for row in itertools.islice(rows, max_samples):
                    w.writerow(row)
                if row is None:
                    return out
            out.append(sio.getvalue())
            sio.close()
