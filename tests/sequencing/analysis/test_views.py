
import pytest

import io
import openpyxl

from tests.sequencing.analysis.adamiya.conftest import *
from tests.sequencing.analysis.adamiya.conftest import _chain_histogram_entry_reads, _chain_amplicon_reads

def test_summary_table(adam_histogram_entry_reads_d, loggedin_client):
    resp = loggedin_client.get("/analysis/TestRun/summary.xlsx")
    assert resp.status_code == 200
    assert resp["Content-Disposition"] == \
        'attachment; filename="TestRun_mapped_reads.xlsx"'
    bio = io.BytesIO()
    for s in resp.streaming_content:
        bio.write(s)
    bio.seek(0)
    # NOTE: filename can also be file-like.
    wb = openpyxl.load_workbook(filename=bio)
    assert wb.active["A1"].value == "Locus"
    # TODO: write the test.
    # for (l_id, bc, inc, amp, msgs), her in adam_histogram_entry_reads_d.items():
