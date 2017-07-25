import pytest
import datetime
import gzip
import os

from django.db.models.query import QuerySet
from sequencing.runs.demux import generate_sample_sheets
from tests.sequencing.analysis.adamiya.conftest import *

@pytest.mark.django_db
def test_get_samplesheet(ngskit, demultiplexingscheme, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a):
    d = datetime.date.today()
    sample_sheet = generate_sample_sheets(
        bcs=[magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a],
        name="Ooga",
        date=d,
        description="Booga",
        kit=ngskit,
        demux_scheme=demultiplexingscheme,
    )
    assert sample_sheet == \
"""[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Experiment Name,Ooga,,,,,,,,
Date,{},,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,Booga,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
151,,,,,,,,,
151,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,1,,,D710,TCCGCGAA,D508,GTACTGAC,,
2,2,,,D718,TGGGAGCC,D502,ATAGAGGC,,
""".format(d)


@pytest.mark.django_db
def test_get_samplesheet_nextseq(ngskit, demultiplexingscheme, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a):
    d = datetime.date.today()
    sample_sheet = generate_sample_sheets(
        bcs=[magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a],
        name="Ooga",
        date=d,
        description="Booga",
        kit=ngskit,
        demux_scheme=demultiplexingscheme,
        rev_left_bc=True,
    )
    assert sample_sheet == \
"""[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Experiment Name,Ooga,,,,,,,,
Date,{},,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,Booga,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
151,,,,,,,,,
151,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,1,,,D710,TCCGCGAA,D508,GTCAGTAC,,
2,2,,,D718,TGGGAGCC,D502,GCCTCTAT,,
""".format(d)


def test_run_demux(ngsrun, demultiplexingscheme, magicalpcr1library, magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a, monkeypatch):
    def mock_run_bcl2fastq(bcl_folder, sample_sheet_path, fastq_folder):
        assert bcl_folder == ngsrun.bcl_directory
        with open(sample_sheet_path) as f:
            sample_sheet = f.read()
        assert sample_sheet == \
"""[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Experiment Name,TestRun,,,,,,,,
Date,{},,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,lib1,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
151,,,,,,,,,
151,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,1,,,D710,TCCGCGAA,D508,GTACTGAC,,
2,2,,,D718,TGGGAGCC,D502,ATAGAGGC,,
""".format(datetime.date.today())
        for num, name in enumerate([
                "Undetermined",
                "1",
                "2",
            ]):
            for read in [1, 2]:
                gzip.open(os.path.join(fastq_folder,
                    "{name}_S{num}_R{read}_001.fastq.gz".format(
                        name=name,
                        num=num,
                        read=read,
                    )), "w").close()
    import sequencing.runs.demux
    monkeypatch.setattr(sequencing.runs.demux, "run_bcl2fastq", mock_run_bcl2fastq)
    srs = list(sequencing.runs.demux.run_demux(ngsrun, demultiplexingscheme))
    demux = srs[0].demux
    fastq_folder = os.path.split(srs[0].fastq1)[0]
    assert demux.ngs_run == ngsrun
    assert demux.demux_scheme == demultiplexingscheme
    for sr, name, num, bc in zip(
        srs,
        ["1", "2"],
        [1, 2],
        [magicalpcr1barcodedcontent, magicalpcr1barcodedcontent_a],
    ):
        assert sr.demux == demux
        assert sr.barcoded_content == bc
        assert sr.library == magicalpcr1library
        # FIXME check with some reads.
        assert sr.num_reads == 0
        for fastq, read in [(sr.fastq1, 1), (sr.fastq2, 2)]:
            assert fastq == os.path.join(fastq_folder,
                "{name}_S{num}_R{read}_001.fastq".format(
                    name=name,
                    num=num,
                    read=read,
            ))
        sr.delete()
    os.rmdir(fastq_folder)

def test_merge_srs(sample_reads_d, mergeddemultiplexing):
    import sequencing.runs.demux
    srs = list(sample_reads_d.values())[:10]
    sr_dict = dict()
    for sr in srs:
        sr_dict.setdefault(sr.barcoded_content, list()).append(sr)

    for bc in sr_dict:
        curr_srs = sr_dict[bc]
        lib = curr_srs[0].library
        merged_sr = sequencing.runs.demux.merge_srs(curr_srs, bc)
        s = 0
        text1 = ""
        text2 = ""
        for sr in curr_srs:
            s += sr.num_reads
            text1 += open(sr.fastq1).read()
            #text1 += "\n" #TODO: Check if necessary
            text2 += open(sr.fastq2).read()
            #text2 += "\n" #TODO: Check if necessary
            assert sr.library == merged_sr["library"]
        res_sr = SampleReads.objects.create(
            demux=mergeddemultiplexing,
            barcoded_content=bc,
            library=merged_sr["library"],
            num_reads=merged_sr["num_reads"],
            fastq1=merged_sr["fastq1"],
            fastq2=merged_sr["fastq2"],
        )

        assert s == res_sr.num_reads
        assert text1 == open(res_sr.fastq1).read()
        assert text2 == open(res_sr.fastq2).read()
        assert lib == res_sr.library
        res_sr.delete()
