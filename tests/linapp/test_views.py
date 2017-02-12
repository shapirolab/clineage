import pytest

import io
import csv
import os

from tests.sequencing.analysis.adamiya.conftest import *

from tests.sequencing.analysis.adamiya.conftest import *
from tests.sequencing.analysis.adamiya.conftest import _chain_histogram_entry_reads, _chain_amplicon_reads
from tests.sequencing.runs.conftest import *
from tests.accounts.conftest import *
from tests.sampling.conftest import *
from tests.sequencing.analysis.full_msv.conftest import *
from tests.wet_storage.conftest import *
from tests.linapp.conftest import *
from django.utils.encoding import smart_text
from sampling.models import FACS

def test_partner_csv_ngs_partner(ngsrun, samplelocation_human_cell_with_se, samplelocation_human_cell_no_se, user, loggedin_client, sample_reads_d):
# def test_partner_csv(ngsrun, human_cell_with_se, user, loggedin_client, sample_reads_d):
    for lib in ngsrun.libraries.all():
        prt=lib.subclass.barcoded_contents[0].amplified_content.cell.individual.partner.get_username()

    assert len(lib.subclass.barcoded_contents)>0
    # bc = lib.subclass.barcoded_contents[0]
    # cell = bc.subclass.amplified_content.cell
    # for cell_cont in cell.amplifiedcontent_set.all():
    #     assert len(cell_cont.physical_locations.exclude(plate__name__contains='AAR'))>0
    resp = loggedin_client.get("/CLineage/csv_view/cell_data/partner_name:{}/ngsrun:{}".format(prt,ngsrun.name))

    assert resp.status_code == 200
    assert resp["Content-Disposition"] == \
           'attachment; filename="TestRun_Shlomo_cell_data.csv"'
    sio = io.StringIO()
    for s in resp.content.decode('utf-8'):
        sio.write(str(s))
    sio.seek(0)
    # NOTE: filename can also be file-like.
    csv_dict = csv.DictReader(sio)
    list_csv= list(csv_dict)
    assert len(list_csv)>0
    for i in list_csv:
        bcid= i['Barcoded Content ID']
        barcoded_cont=lib.subclass.barcoded_contents.get(id=bcid)
        cell = barcoded_cont.subclass.amplified_content.cell
        sr = barcoded_cont.samplereads_set.get(library=lib)
        for cell_cont in cell.amplifiedcontent_set.all():
            sampling = cell.sampling
            if sampling:
                try:
                    facs = sampling.facs
                except FACS.DoesNotExist:
                    facs = None
            else:
                facs = None
            for loc in cell_cont.physical_locations.exclude(plate__name__contains='AAR'):
                assert i['Sample Reads ID']== '{}'.format(sr.id)
                assert i['CellContent ID']== smart_text(cell_cont.pk)
                assert i['Cell ID']== smart_text(cell.pk)
                assert i['Cell Name']== smart_text(cell.name)
                assert i['Cell Group']== smart_text(cell.classification)
                assert i['Individual Name']== smart_text(cell.individual.name)
                assert i['Individual Comment']== smart_text(cell.individual.comment)
                assert i['Extraction Event']== smart_text(cell.sampling.extraction.extraction_event.name if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else '')
                assert i['Extraction Event Comment']== smart_text(cell.sampling.extraction.extraction_event.comment if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else '')
                assert i['Gender']== smart_text(cell.individual.sex)
                assert i['Sample Name']== smart_text(cell.sampling.extraction.name if cell.sampling else '')
                assert i['Sample Comment']== smart_text(cell.sampling.extraction.comment if cell.sampling else '')
                assert i['Organ']== smart_text(cell.sampling.extraction.organ.name if cell.sampling else '')
                assert i['Tissue']== smart_text(cell.sampling.extraction.tissue.name if cell.sampling else '')
                assert i['Sampling Event']== smart_text(cell.sampling.name if cell.sampling else '')
                assert i['Sampling Comment']== smart_text(cell.sampling.comment if cell.sampling else '')
                assert i['FACS Marker']== smart_text(facs.marker.name if facs else '')
                assert i['Cell Type']== smart_text(cell.composition.name)
                assert i['Plate']== smart_text(loc.plate.name)
                assert i['Well']== smart_text(loc.well)
                assert i['Plate Location']== smart_text(loc.plate.platestorage_set.all()[0] if loc.plate.platestorage_set.all() else '')


                #assert cct["A1"].value == "Barcoded_Content_ID"
    # TODO: write the test.
    # for (l_id, bc, inc, amp, msgs), her in adam_histogram_entry_reads_d.items():

def test_partner_csv(ngsrun, samplelocation_human_cell_with_se, samplelocation_human_cell_no_se, user, loggedin_client, sample_reads_d):
# def test_partner_csv(ngsrun, human_cell_with_se, user, loggedin_client, sample_reads_d):
    for lib in ngsrun.libraries.all():
        prt=lib.subclass.barcoded_contents[0].amplified_content.cell.individual.partner.get_username()
        ind = lib.subclass.barcoded_contents[0].amplified_content.cell.individual.name
    assert len(lib.subclass.barcoded_contents)>0
    # bc = lib.subclass.barcoded_contents[0]
    # cell = bc.subclass.amplified_content.cell
    # for cell_cont in cell.amplifiedcontent_set.all():
    #     assert len(cell_cont.physical_locations.exclude(plate__name__contains='AAR'))>0
    resp = loggedin_client.get("/CLineage/csv_view/cell_data/")
    inside_content=()
    field_names=("/partner_name:{}".format(prt), "/individual_name:{}".format(ind), "/ngsrun:{}".format(ngsrun.name))
    url_positions=[(field_names[0],field_names[1],field_names[2]),(field_names[0],field_names[1],''),(field_names[0],'',field_names[2]),(field_names[0],'',''),('', field_names[1],field_names[2]),('', field_names[1], ''),( '', '' , field_names[2]),('','','')]
    for partner_option, individual_option, ngs_option in url_positions:
        resp = loggedin_client.get("/CLineage/csv_view/cell_data{}{}{}".format(partner_option,individual_option,ngs_option))

        assert resp.status_code == 200
        # assert resp["Content-Disposition"] == \
        #        'attachment; filename="TestRun_Shlomo_cell_data.csv"'
        sio = io.StringIO()
        for s in resp.content.decode('utf-8'):
            sio.write(str(s))
        sio.seek(0)
        # NOTE: filename can also be file-like.
        csv_dict = csv.DictReader(sio)
        list_csv= list(csv_dict)
        assert len(list_csv)>0
        for i in list_csv:
            bcid= i['Barcoded Content ID']
            barcoded_cont=lib.subclass.barcoded_contents.get(id=bcid)
            cell = barcoded_cont.subclass.amplified_content.cell
            sr = barcoded_cont.samplereads_set.get(library=lib)
            for cell_cont in cell.amplifiedcontent_set.all():
                sampling = cell.sampling
                if sampling:
                    try:
                        facs = sampling.facs
                    except FACS.DoesNotExist:
                        facs = None
                else:
                    facs = None
                for loc in cell_cont.physical_locations.exclude(plate__name__contains='AAR'):
                    assert i['Sample Reads ID']== '{}'.format(sr.id)
                    assert i['CellContent ID']== smart_text(cell_cont.pk)
                    assert i['Cell ID']== smart_text(cell.pk)
                    assert i['Cell Name']== smart_text(cell.name)
                    assert i['Cell Group']== smart_text(cell.classification)
                    assert i['Individual Name']== smart_text(cell.individual.name)
                    assert i['Individual Comment']== smart_text(cell.individual.comment)
                    assert i['Extraction Event']== smart_text(cell.sampling.extraction.extraction_event.name if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else '')
                    assert i['Extraction Event Comment']== smart_text(cell.sampling.extraction.extraction_event.comment if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else '')
                    assert i['Gender']== smart_text(cell.individual.sex)
                    assert i['Sample Name']== smart_text(cell.sampling.extraction.name if cell.sampling else '')
                    assert i['Sample Comment']== smart_text(cell.sampling.extraction.comment if cell.sampling else '')
                    assert i['Organ']== smart_text(cell.sampling.extraction.organ.name if cell.sampling else '')
                    assert i['Tissue']== smart_text(cell.sampling.extraction.tissue.name if cell.sampling else '')
                    assert i['Sampling Event']== smart_text(cell.sampling.name if cell.sampling else '')
                    assert i['Sampling Comment']== smart_text(cell.sampling.comment if cell.sampling else '')
                    assert i['FACS Marker']== smart_text(facs.marker.name if facs else '')
                    assert i['Cell Type']== smart_text(cell.composition.name)
                    assert i['Plate']== smart_text(loc.plate.name)
                    assert i['Well']== smart_text(loc.well)
                    assert i['Plate Location']== smart_text(loc.plate.platestorage_set.all()[0] if loc.plate.platestorage_set.all() else '')


                    #assert cct["A1"].value == "Barcoded_Content_ID"
        # TODO: write the test.
        # for (l_id, bc, inc, amp, msgs), her in adam_histogram_entry_reads_d.items():
