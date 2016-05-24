import pytest
import datetime

from django.db.models.query import QuerySet
from sequencing.runs.demux import generate_sample_sheets

@pytest.mark.django_db
def test_get_samplesheet(ngsrun,demultiplexingscheme):
    # TODO: sort out adaptors
    bc_l = []
    for library in ngsrun.libraries.select_subclasses():
        it = library.barcoded_contents
        if isinstance(it, QuerySet):
            it = it.select_related(
                "barcodes",
                "barcodes__left",
                "barcodes__right"
            )
        for bc in it:
            bc_l.append(bc)
    description = "_".join(lib.name for lib in ngsrun.libraries.all())
    sample_sheet = generate_sample_sheets(
        bcs=bc_l,
        name=ngsrun.name,
        date=ngsrun.date,
        description=description,
        kit=ngsrun.kit,
        demux_scheme=demultiplexingscheme,
    )
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

