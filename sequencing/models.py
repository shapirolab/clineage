from django.db import models

from lib_prep.workflows.models import CellContent
from sequencing.runs.models import MergedReads
from targeted_enrichment.planning.models import Target


### -------------------------------------------------------------------------------------
class SequencingData(models.Model): # This contains the actual data.
    cell_content = models.ForeignKey(CellContent)
    merged_reads = models.ForeignKey(MergedReads)
    target = models.ForeignKey(Target)
    target_offset = models.IntegerField(null=True)
    fastq = models.FilePathField(null=True)
    vcf = models.FilePathField(null=True)

    class Meta:
        index_together = (
            ("cell_content", "merged_reads", "target")
        )
