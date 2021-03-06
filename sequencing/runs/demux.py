
from plumbum import local
import csv
import io
import os
import resource
import shutil
from Bio import SeqIO

from misc.utils import get_unique_path, unique_dir_cm, unlink, unique_file_cm

from django.db.models.query import QuerySet

from sequencing.runs.models import Demultiplexing, MergedDemultiplexing
from sequencing.analysis.models import SampleReads
from plumbum import ProcessExecutionError


def run_bcl2fastq(bcl_folder, sample_sheet_path, fastq_folder):
    assert resource.getrlimit(resource.RLIMIT_OFILE)[0] >= 4096
    assert resource.getrlimit(resource.RLIMIT_OFILE)[1] >= 4096
    # resource.setrlimit(resource.RLIMIT_OFILE, (4096, 4096))
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


def generate_sample_sheets(bcs, name, date, description, kit, demux_scheme, rev_left_bc=False):
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
            LEFT_BARCODE_SEQ: bc.barcodes.left.ref_sequence if not rev_left_bc else bc.barcodes.left.ref_sequence.rev_comp(),
            RIGHT_BARCODE_ID: bc.barcodes.right.name,
            RIGHT_BARCODE_SEQ: bc.barcodes.right.ref_sequence,
        })
    s = sio.getvalue()
    sio.close()
    return s


UNDETERMINED = "Undetermined"
SAMPLE_FASTQ_GZ_FORMAT = "{id}_S{idx}_R{read:d}_001.fastq.gz"
SAMPLE_FASTQ_FORMAT = "{id}_S{idx}_R{read:d}_001.fastq"


def run_demux(ngs_run, demux_scheme, write_her_files=False):
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
        rev_left_bc=ngs_run.machine.type.rev_left_bc,
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
                try:
                    local["gunzip"](fastqgz)
                except ProcessExecutionError as e:
                    # TODO: assert 0 reads against demux reports
                    if e.args[3].find('No such file or directory'):  # demux didn't output the file due to no reads
                        for read in [1, 2]:  # Manually touch the files for downstream analysis
                            with open(os.path.join(fastq_folder,
                                      SAMPLE_FASTQ_FORMAT.format(
                                          id=bc.id,
                                          idx=idx,
                                          read=read,
                                      )), 'w'):
                                pass  # Nothing to do here but touch the file
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
                write_her_files=write_her_files,
            )


def merge_srs(srs_lst, bc, root_demux):
    """ This method gets list of SampleReads from same Barcoded Content and merges them into one Sample Read"""
    try:
        yield SampleReads.objects.get(
            demux=root_demux,
            barcoded_content=bc)
    except SampleReads.DoesNotExist:
        merged_sr = dict()
        merged_sr["num_reads"] = 0
        with unique_file_cm("fastq") as fastq1_merged:
            with open(fastq1_merged, 'w') as fd1:
                with unique_file_cm("fastq") as fastq2_merged:
                    with open(fastq2_merged, 'w') as fd2:
                        keep_her_file = False
                        for sr in srs_lst:
                            keep_her_file = keep_her_file or sr.write_her_files # if one keep hers so all will write
                            assert sr.barcoded_content == bc
                            merged_sr["num_reads"] += sr.num_reads
                            with open(sr.fastq1, 'r') as fastq1_in:
                                shutil.copyfileobj(fastq1_in, fd1)
                            with open(sr.fastq2, 'r') as fastq2_in:
                                shutil.copyfileobj(fastq2_in, fd2)
                        merged_sr["fastq1"] = fastq1_merged
                        merged_sr["fastq2"] = fastq2_merged
                        merged_sr["library"] = sr.library #bc-library connections is 1 to 1 so this assignment is OK
                        merged_sr["bc"] = bc
                        merged_sr["write_her_files"] = keep_her_file
        yield SampleReads.objects.create(
            demux=root_demux,
            barcoded_content=bc,
            library=merged_sr["library"],
            num_reads=merged_sr["num_reads"],
            fastq1=merged_sr["fastq1"],
            fastq2=merged_sr["fastq2"],
            write_her_files=merged_sr['write_her_files']
        )


def prepare_merge_demuxes(ngs_runs, merged_demux_scheme):
    """
    This method gets list of runs and creates dictionary with sample reads to merge according to barcoded content
    Each run must have its demux ready!!
    """

    assert len(ngs_runs) > 1  # cannot merge one run

    root_demux = MergedDemultiplexing.objects.create(
        ngs_run=ngs_runs[0],
        demux_scheme=merged_demux_scheme,
    )
    # workaround to bypass voodoo of unknown field
    root_demux.ngs_runs = ngs_runs[1:]
    root_demux.save()

    srs_by_bc_dict = dict()
    # populate a dictionary with barcoded content id as key and list of its Sample Reads as value
    for ngs_run in ngs_runs:
        demux = [dm for dm in ngs_run.demultiplexing_set.all() if not MergedDemultiplexing.objects.filter(id=dm.id)]
        assert len(demux) == 1
        demux = demux[0]
        for sr in demux.samplereads_set.all():
            srs_by_bc_dict.setdefault(sr.barcoded_content, list()).append(sr)

    if merged_demux_scheme.name == 'intersection_demux_scheme':  # filter out BCs that are only present in a single run
        return root_demux, {bc: srs for bc, srs in srs_by_bc_dict.items() if len(srs) == len(ngs_runs)}
    return root_demux, srs_by_bc_dict


def merge_demuxes(ngs_runs, merged_demux_scheme):
    """
    This method gets list of runs and merges its demultiplexing.
    Meaning, it concatenates the fastq files into one
    Each run must have its demux ready!!
    """
    root_demux, srs_by_bc_dict = prepare_merge_demuxes(ngs_runs, merged_demux_scheme)

    # create merged file for each fastq1,fastq2 of each Sample Read
    for bc in srs_by_bc_dict:
        srs = srs_by_bc_dict[bc]
        merged_sr = merge_srs(srs, bc, root_demux)
        yield merged_sr


def merge_demuxes_parallel(executor, root_demux, srs_by_bc_dict):
    """
    Like merge_demuxes but in parallel
    Each run must have its demux ready!!
    """
    for bc in srs_by_bc_dict:
        srs = srs_by_bc_dict[bc]
        merged_sr = executor.submit(merge_srs, srs, bc, root_demux)
        yield merged_sr
