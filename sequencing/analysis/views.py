
from io import BytesIO
import xlsxwriter
from Bio import SeqIO

from django.http import FileResponse, HttpResponseNotFound
from django.db.models import Sum

from sequencing.runs.models import NGSRun, Demultiplexing
from sequencing.analysis.models import SampleReads, Histogram, \
    HistogramEntryReads

def summary_table(request, ngs_run):
    try:
        nr = NGSRun.objects.get(name=ngs_run)
    except NGSRun.DoesNotExist:
        return HttpResponseNotFound(b"Sorry, run does not exist.")
    try:
        # FIXME
        demux = Demultiplexing.objects.get(ngs_run__name=ngs_run)
    except Demultiplexing.DoesNotExist:
        return HttpResponseNotFound(b"Sorry, run not analysed.")
    bio = BytesIO()
    workbook = xlsxwriter.Workbook(bio)
    for library in nr.libraries.select_subclasses():
        bcs = list(library.barcoded_contents)
        # rev_bcs = {bc: i+1 for i in enumerate(bcs)}
        srs = list(SampleReads.objects.filter(
            demux=demux,
            library=library,
            barcoded_content__in=bcs,
        ))
        assert len(srs) == len(bcs)
        rev_sr_ids = {sr.id: i+1 for i, sr in enumerate(srs)}
        amplicons = list(library.amplicons)
        rev_amp_ids = {amp.id: i+1 for i, amp in enumerate(amplicons)}
        SUMMARY_ROWS = 1 + len(amplicons)

        w = workbook.add_worksheet(name=library.name)

        w.write(0, 0, "Locus")
        w.write(SUMMARY_ROWS, 0, "TotalReads")
        for sr in srs:
            col = rev_sr_ids[sr.id]
            w.write(0, col, sr.barcoded_content.subclass.amplified_content.name)
            # FIXME
            with open(sr.fastq1) as f:
                seq = SeqIO.parse(f, format="fastq")
                n = sum(1 for x in seq)
            w.write(SUMMARY_ROWS, col, n)
        for amp_id, row in rev_amp_ids.items():
            # TODO: is this the id we want? teid? ter-id?
            w.write(row, 0, amp_id)

        for h in Histogram.objects.filter(
            # FIXME: microsatellites_version
            amplicon__in=amplicons,
            sample_reads__in=srs,
        ):
            n = HistogramEntryReads.objects.filter(histogram=h) \
                .aggregate(Sum('num_reads')).popitem()[1]
            w.write(
                rev_amp_ids[h.amplicon_id],
                rev_sr_ids[h.sample_reads_id],
                n,
            )
    workbook.close()
    bio.seek(0)
    response = FileResponse(bio,
        content_type="application/vnd.openxmlformats-officedocument" \
            ".spreadsheetml.sheet"
        )
    response["Content-Disposition"] = \
        'attachment; filename="{}_mapped_reads.xlsx"'.format(nr.name)
    return response
