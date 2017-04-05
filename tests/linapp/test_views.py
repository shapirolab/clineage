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

from tests.targeted_enrichment.planning.conftest import *
from tests.lib_prep.multiplexes.conftest import *


from django.utils.encoding import smart_text
from sampling.models import FACS
from lib_prep.multiplexes.models import *
from targeted_enrichment.reagents.models import *


def test_partner_csv_no_physical_location(ngsrun, amplifiedcontent, user, loggedin_client, sample_reads_d):
    for lib in ngsrun.libraries.all():
        prt=lib.subclass.barcoded_contents[0].amplified_content.cell.individual.partner.get_username()
        ind = lib.subclass.barcoded_contents[0].amplified_content.cell.individual.name
        assert len(lib.subclass.barcoded_contents)>0
        resp = loggedin_client.get("/CLineage/csv_view/cell_data/")
        field_names=("/partner_name:{}".format(prt), "/individual_name:{}".format(ind), "/ngsrun:{}".format(ngsrun.name))
        url_positions=[(field_names[0],field_names[1],field_names[2]),(field_names[0],field_names[1],''),(field_names[0],'',field_names[2]),(field_names[0],'',''),('', field_names[1],field_names[2]),('', field_names[1], ''),( '', '' , field_names[2]),('','','')]
        for partner_option, individual_option, ngs_option in url_positions:
            resp = loggedin_client.get("/CLineage/csv_view/cell_data{}{}{}".format(partner_option,individual_option,ngs_option))
            assert resp.status_code == 200
            sio = io.StringIO()
            for s in resp.content.decode('utf-8'):
                sio.write(str(s))
            sio.seek(0)
            # NOTE: filename can also be file-like.
            csv_dict = csv.DictReader(sio)
            list_csv= list(csv_dict)
            assert len(list_csv)>0
            assert amplifiedcontent.id in set(int(i['Barcoded Content ID']) for i in list_csv)


def test_partner_csv(ngsrun, samplelocation_human_cell_with_se, samplelocation_human_cell_no_se, user,
                     loggedin_client, sample_reads_d):
    for lib in ngsrun.libraries.all():
        prt = lib.subclass.barcoded_contents[0].amplified_content.cell.individual.partner.get_username()
        ind = lib.subclass.barcoded_contents[0].amplified_content.cell.individual.name
        assert len(lib.subclass.barcoded_contents) > 0
        resp = loggedin_client.get("/CLineage/csv_view/cell_data/")
        field_names = (
        "/partner_name:{}".format(prt), "/individual_name:{}".format(ind), "/ngsrun:{}".format(ngsrun.name))
        url_positions = [(field_names[0], field_names[1], field_names[2]),
                         (field_names[0], field_names[1], ''), (field_names[0], '', field_names[2]),
                         (field_names[0], '', ''), ('', field_names[1], field_names[2]),
                         ('', field_names[1], ''), ('', '', field_names[2]), ('', '', '')]
        for partner_option, individual_option, ngs_option in url_positions:
            resp = loggedin_client.get(
                "/CLineage/csv_view/cell_data{}{}{}".format(partner_option, individual_option, ngs_option))
            assert resp.status_code == 200
            sio = io.StringIO()
            for s in resp.content.decode('utf-8'):
                sio.write(str(s))
            sio.seek(0)
            # NOTE: filename can also be file-like.
            csv_dict = csv.DictReader(sio)
            list_csv = list(csv_dict)
            assert len(list_csv) > 0
            for i in list_csv:
                bcid = i['Barcoded Content ID']
                barcoded_cont = lib.subclass.barcoded_contents.get(id=bcid)
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
                        assert i['Sample Reads ID'] == '{}'.format(sr.id)
                        assert i['CellContent ID'] == smart_text(cell_cont.pk)
                        assert i['Cell ID'] == smart_text(cell.pk)
                        assert i['Cell Name'] == smart_text(cell.name)
                        assert i['Cell Group'] == smart_text(cell.classification)
                        assert i['Individual Name'] == smart_text(cell.individual.name)
                        assert i['Individual Comment'] == smart_text(cell.individual.comment)
                        assert i['Extraction Event'] == smart_text(
                            cell.sampling.extraction.extraction_event.name if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else '')
                        assert i['Extraction Event Comment'] == smart_text(
                            cell.sampling.extraction.extraction_event.comment if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else '')
                        assert i['Gender'] == smart_text(cell.individual.sex)
                        assert i['Sample Name'] == smart_text(
                            cell.sampling.extraction.name if cell.sampling else '')
                        assert i['Sample Comment'] == smart_text(
                            cell.sampling.extraction.comment if cell.sampling else '')
                        assert i['Organ'] == smart_text(
                            cell.sampling.extraction.organ.name if cell.sampling else '')
                        assert i['Tissue'] == smart_text(
                            cell.sampling.extraction.tissue.name if cell.sampling else '')
                        assert i['Sampling Event'] == smart_text(
                            cell.sampling.name if cell.sampling else '')
                        assert i['Sampling Comment'] == smart_text(
                            cell.sampling.comment if cell.sampling else '')
                        assert i['FACS Marker'] == smart_text(facs.marker.name if facs else '')
                        assert i['Cell Type'] == smart_text(cell.composition.name)
                        assert i['Plate'] == smart_text(loc.plate.name)
                        assert i['Well'] == smart_text(loc.well)
                        assert i['Plate Location'] == smart_text(
                            loc.plate.platestorage_set.all()[0] if loc.plate.platestorage_set.all() else '')


#padlockpanel
def test_panel_csv(user, loggedin_client, pcr1panel, padlockpanel, requires_microsatellites, requires_snps):
    panels=[pcr1panel, padlockpanel]
    for panel in panels:
        panel_name=panel.name if ' ' not in panel.name else panel.name.split(' ')[1]
        field_names = "/panel:{}".format(panel_name)
        for panel_option in [field_names, '/','']:
            resp = loggedin_client.get("/CLineage/csv_view/panel_data{}".format(panel_option))
            assert resp.status_code == 200
            sio = io.StringIO()
            for s in resp.content.decode('utf-8'):
                sio.write(str(s))
            sio.seek(0)
            # NOTE: filename can also be file-like.
            csv_dict = csv.DictReader(sio)
            list_csv = list(csv_dict)
            assert len(list_csv) > 0
            for i in list_csv:
                terid = i['TER ID']
                if i['Mpx groups names' ] is not '':
                    ter = PCR1PrimerPairTER.objects.get(id=terid)
                else:
                    ter = OM6PadlockTER.objects.get(id=terid)
                amplicon = ter.amplicon
                target_id = i['Target ID']
                snp=None
                ms=None
                for slice_contains in amplicon.slice.contains.all():
                    if slice_contains.target_set.filter(id=target_id):
                        target = slice_contains.target_set.get(id=target_id)
                try:
                    ms = target.microsatellite
                except Microsatellite.DoesNotExist:
                    try:
                        snp = target.snp
                    except SNP.DoesNotExist:
                        assert 0==1, 'not a snp nor a ms???'

                assert i['MS ID'] == smart_text('{}'.format(ms.id) if ms is not None else '')
                assert i['SNP ID'] == smart_text('{}'.format(snp.id) if snp is not None else '')
                assert i ['Target ID' ] == smart_text(target.id if target is not None else '')
                assert i ['Old Adam TEpk' ] == smart_text(ter.old_adam_te_pk)
                assert i ['TER ID' ] == smart_text(ter.id)
                assert i ['Amplicon ID' ] == smart_text(amplicon.id)
                assert i ['Basic Unit size' ] == smart_text(ms.repeat_unit_len if ms is not None else '')
                assert i ['Expected Number of repeats' ] == smart_text(ms.repeat_number if ms is not None else '')
                assert i ['Basic Unit Type' ] == smart_text(ms.repeat_unit_type if ms is not None else '')
                assert i ['SNP Sequence' ] ==smart_text(snp.slice.sequence.decode('utf-8') if snp is not None else '')
                assert i ['SNP modified' ] ==smart_text(snp.modified if snp is not None else '')
                assert i ['Chromosome' ] == smart_text(ms.slice.chromosome if ms is not None else snp.slice.chromosome)
                assert i ['Length MS' ] == smart_text(abs(ms.slice.end_pos - ms.slice.start_pos) if ms is not None else '')
                assert i ['Target sequence' ] == smart_text(target.slice.sequence.seq.decode('utf-8'))
                assert i ['Primer sequence - Fw' ] == smart_text(ter.left_primer.sequence.seq.decode('utf-8') \
                                                                        if i[
                        'Mpx groups names'] is not '' else ter.padlock.left_ugs.sequence.seq.decode('utf-8'))
                assert i ['Primer sequence -  Rev' ] == smart_text(ter.right_primer.sequence.seq.decode('utf-8') \
                                                                        if i[
                        'Mpx groups names'] is not '' else ter.padlock.right_ugs.sequence.seq.decode('utf-8'))
                assert i ['Target location on Chromosome - start' ] == smart_text(ter.te.left.slice.start_pos)
                assert i ['Target location on Chromosome - end' ] == smart_text(ter.te.right.slice.end_pos)
                assert i ['Fw primer location on Chromosome - start' ] == smart_text(ter.left_primer.ugs.slice.start_pos \
                                                                        if i[
                        'Mpx groups names'] is not '' else ter.padlock.left_ugs.slice.start_pos)
                assert i ['Fw primer location on Chromosome - end' ] == smart_text(ter.left_primer.ugs.slice.end_pos \
                                                                        if i[
                        'Mpx groups names'] is not '' else ter.padlock.left_ugs.slice.end_pos)
                assert i ['Rev primer location on Chromosome - start' ] == smart_text(ter.right_primer.ugs.slice.start_pos \
                                                                        if i[
                        'Mpx groups names'] is not '' else ter.padlock.right_ugs.slice.start_pos)
                assert i ['Rev primer location on Chromosome - end' ] == smart_text(ter.right_primer.ugs.slice.end_pos \
                                                                        if i[
                        'Mpx groups names'] is not '' else ter.padlock.right_ugs.slice.end_pos)
                assert i ['Amplicon location on Chromosome - start' ] == smart_text(amplicon.slice.start_pos, )
                assert i ['Amplicon location on Chromosome - end' ] == smart_text(amplicon.slice.end_pos)
                assert i ['Amplicon Sequence' ] == smart_text(amplicon.slice.sequence.seq.decode('utf-8'))
                assert i ['MS location on Chromosome - start' ] == smart_text(ms.slice.start_pos if ms is not None else '')
                assert i ['MS location on Chromosome - end' ] == smart_text(ms.slice.end_pos if ms is not None else '')
                assert i ['SNP location on the Chromosome' ] ==smart_text(snp.slice.start_pos if snp is not None else '')
                # assert i ['Mpx groups names' ] == smart_text(mx.name if type(panel) == PCR1Panel else '')
                # assert i ['Mpx group size' ] == smart_text(len(mx.ters.select_subclasses()) if type(panel) == PCR1Panel else '')
                # assert i ['Mixs groups names' ] ==smart_text(mx.name if type(panel) == OM6Panel else '')
                # assert i ['Mixs group size' ] == smart_text(len(mx.ters.select_subclasses()) if type(panel) == OM6Panel else '')
