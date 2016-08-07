import os
import csv
import seaborn as sns  # development
from collections import defaultdict

from django.utils.encoding import smart_text
from django.contrib.auth.models import User

from linapp.models import UserReport
from sampling.models import Cell, Individual, FACS


def query_partner_individuals(partner_name, individual_name=None):
    partner = User.objects.get(username__contains=partner_name)
    if individual_name:
        individuals = Individual.objects.filter(partner=partner).filter(name__contains=individual_name)
    else:
        individuals = Individual.objects.filter(partner=partner)
    return partner, individuals


def to_hex(n):
    return hex(int(n*255))[2:].upper()


def hex_to_rgb(color_map, cell):
    rgb = []
    color_cell = color_map[cell]
    for i in color_cell:
        rgb.append(to_hex(i))
    return rgb


def extraction_events_for_individual(indi):
    for ee in sorted(list(indi.extractionevent_set.all()), key=lambda i: i.name):
        yield ee


def extractions_for_individual(indi):
    for ee in extraction_events_for_individual(indi):
        for e in sorted(list(ee.extraction_set.all()), key=lambda i: i.name):
            yield e, ee


def sampling_event_for_individual(indi, report_dict):
    for e, ee in extractions_for_individual(indi):
        report_dict[indi.name][ee.name][e.name]['Extraction_date'] = ee.date
        if e.samplingevent_set.all():
            for se in sorted(list(e.samplingevent_set.all()), key=lambda i: i.name):
                yield se, e, ee
        else:
            yield None, e, ee


def sorted_cell_classifications(cells):
    return sorted(list(set([cell.classification for cell in cells])), key=lambda x: (x is None, x))


def cells_for_individual(indi, report_dict):
    for se, e, ee in sampling_event_for_individual(indi, report_dict):
        if se:
            for cls in sorted_cell_classifications(se.cell_set.all()):
                if se.cell_set.filter(classification=cls):
                    yield cls, se, e, ee


def populate_report_dict(indi, report_dict, cellrow, color_map):
    for cls, se, e, ee in cells_for_individual(indi, report_dict):
        report_dict[indi.name][ee.name][e.name][se.name][str(cls)]['Cells_separation_date'] = se.date
        report_dict[indi.name][ee.name][e.name][se.name][str(cls)]['Cells_separation_details'] = se.comment
        report_dict[indi.name][ee.name][e.name][se.name][str(cls)]['Cells_classification_string'] = cls
        report_dict[indi.name][ee.name][e.name][se.name][str(cls)]['Cells_color'] = hex_to_rgb(color_map, se.cell_set.filter(classification=cls)[0])
        report_dict[indi.name][ee.name][e.name][se.name][str(cls)]['Cells_pos'] = [cellrow+1, cellrow+se.cell_set.filter(classification=cls).count()]
        cellrow += se.cell_set.filter(classification=cls).count()
    return report_dict, cellrow


def collect_cells_without_clasification(report_dict, indi, cellrow, color_map):
    report_dict[indi.name]['cells_list'] = 1 #temporary workaround TODO:revise
    for cls in sorted_cell_classifications(indi.cell_set.all()):
        report_dict[indi.name][str(cls)]['Cells_color'] = hex_to_rgb(color_map, indi.cell_set.filter(classification=cls)[0])
        report_dict[indi.name][str(cls)]['Cells_pos'] = [cellrow+1, cellrow+indi.cell_set.filter(classification=cls).count()]
        report_dict[indi.name][str(cls)]['Cells_classification_string'] = cls
        cellrow += indi.cell_set.filter(classification=cls).count()
    return report_dict, cellrow

def get_partner_report(partner_name, individual_name=None, palette_name='hls'):
    cellrow = 1
    report_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))))
    color_map = get_cells_color_map(get_cells_grouping(partner_name, individual_name), palette_name)
    partner, individuals = query_partner_individuals(partner_name, individual_name)
    report_dict['Collaborator'] = partner.username
    cells = get_cells(partner_name, individual_name)
    ur = UserReport.get_create_new(cells=cells, partner=partner, individual=individuals)
    report_dict['ID'] = ur.pk
    for individual in sorted(list(individuals), key=lambda i: i.name):
        report_dict[partner.username][individual.name]['name'] = individual.name
        if not individual.extractionevent_set.all():
            report_dict[partner.username], cellrow = collect_cells_without_clasification(report_dict[partner.username], individual, cellrow, color_map)
        report_dict[partner.username][individual.name]['Collaborator_table'] = None
        if individuals.count() > 1:
            report_dict[partner.username][individual.name]['Database_table'] = '{}_cell_data.csv'.format(partner.username)
        else:
            report_dict[partner.username][individual.name]['Database_table'] = '{}_{}_cell_data.csv'.format(partner.username, individual.name)
        report_dict[partner.username], cellrow = populate_report_dict(individual, report_dict[partner.username], cellrow, color_map)
    return report_dict


def get_cell_files_from_folder(cell_folder):
    folder_cells = {}
    for file_name in os.listdir(cell_folder):
        if ('hist_' in file_name) or (not 'hist' in file_name):
            continue
        folder_cells[Cell.objects.get(id=int(file_name.split('-')[0]))] = file_name
    return folder_cells


def get_cells_filenames_in_folder(cells, cell_folder):
    folder_cells = get_cell_files_from_folder(cell_folder) if cell_folder else None
    for cell in cells:
        if folder_cells and cell in folder_cells:
            yield cell, folder_cells[cell]
        else:
            yield cell, None


def get_cells(partner_name, individual_name=None, cell_groups=None):
    cells = set()
    if not cell_groups:
        cell_groups = get_cells_grouping(partner_name, individual_name=individual_name)
    cells |= set(cell_groups.keys())
    partner, individuals = query_partner_individuals(partner_name, individual_name=None)
    for i in individuals:
        cells |= set(i.cell_set.all())
    return list(cells)


#TODO: DRY
def sorted_cells(individual):
    cell_list = []
    if not individual.extractionevent_set.all() or \
            not individual.extractionevent_set.all()[0].extraction_set.all() or \
            not individual.extractionevent_set.all()[0].extraction_set.all()[0].samplingevent_set.all():
        for cls in sorted_cell_classifications(individual.cell_set.all()):
            for cell in individual.cell_set.filter(classification=cls):
                cell_list.append(cell)
    for ee in sorted(list(individual.extractionevent_set.all()), key=lambda ee: ee.name):
        for e in sorted(list(ee.extraction_set.all()), key=lambda e: e.name):
            for se in sorted(list(e.samplingevent_set.all()), key=lambda se: se.name):
                for cls in sorted_cell_classifications(se.cell_set.all()):
                    for cell in se.cell_set.filter(classification=cls):
                        cell_list.append(cell)
    return cell_list


def get_cells_grouping(partner_name, individual_name=None, current_group=0):
    partner, individuals = query_partner_individuals(partner_name, individual_name)
    cell_groups = {}
    for individual in sorted(list(individuals), key=lambda i: i.name):
        if not individual.extractionevent_set.all() or \
                not individual.extractionevent_set.all()[0].extraction_set.all() or \
                not individual.extractionevent_set.all()[0].extraction_set.all()[0].samplingevent_set.all():
            for cls in sorted_cell_classifications(individual.cell_set.all()):
                for cell in individual.cell_set.filter(classification=cls):
                    cell_groups[cell] = current_group
                if individual.cell_set.filter(classification=cls):
                    current_group += 1
            continue
        for ee in sorted(list(individual.extractionevent_set.all()), key=lambda ee: ee.name):
            for e in sorted(list(ee.extraction_set.all()), key=lambda e: e.name):
                for se in sorted(list(e.samplingevent_set.all()), key=lambda se: se.name):
                    for cls in sorted_cell_classifications(se.cell_set.all()):
                        for cell in se.cell_set.filter(classification=cls):
                            cell_groups[cell] = current_group
                        if se.cell_set.filter(classification=cls):
                            current_group += 1
    return cell_groups


def get_cells_color_map(cell_groups, palette_name='hls'):
    if len(set(cell_groups.values())) == 1:
        return {cell: (0, 0.5, 0.5) for cell in cell_groups}
    palette = sns.color_palette(palette_name, len(set(cell_groups.values())))  # development
    return {cell: palette[cell_groups[cell]] for cell in cell_groups}  # development
    # return {cell: (0, 0.5, 0.5) for cell in cell_groups}


def user_cells_table_values(partner_name, individual_name=None, cell_folder=None, palette_name='hls'):
    partner, individuals = query_partner_individuals(partner_name, individual_name)
    color_map = get_cells_color_map(get_cells_grouping(partner_name, individual_name), palette_name)
    for individual in sorted(list(individuals), key=lambda i: i.name):
        for cell, file_name in get_cells_filenames_in_folder(sorted_cells(individual), cell_folder):
            for cell_cont in cell.amplifiedcontent_set.all():
                # assert cell_cont.physical_locations.exclude(plate__name__contains='AAR').count() <= 1
                sampling = cell.sampling
                if sampling:
                    try:
                        facs = sampling.facs
                    except FACS.DoesNotExist:
                        facs = None
                else:
                    facs = None
                for loc in cell_cont.physical_locations.exclude(plate__name__contains='AAR'):
                    yield {
                        'CellContent ID': smart_text(cell_cont.pk),
                        'Cell ID': smart_text(cell.pk),
                        'Sequencing File Name': smart_text(file_name),
                        'Cell Name': smart_text(cell.name),
                        'Cell Group': smart_text(cell.classification),
                        'Individual Name': smart_text(cell.individual.name),
                        'Individual Comment': smart_text(cell.individual.comment),
                        'Extraction Event': smart_text(cell.sampling.extraction.extraction_event.name if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else ''),
                        'Extraction Event Comment': smart_text(cell.sampling.extraction.extraction_event.comment if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else ''),
                        'Gender': smart_text(cell.individual.sex),
                        'Sample Name': smart_text(cell.sampling.extraction.name if cell.sampling else ''),
                        'Sample Comment': smart_text(cell.sampling.extraction.comment if cell.sampling else ''),
                        'Organ': smart_text(cell.sampling.extraction.organ.name if cell.sampling else ''),
                        'Tissue': smart_text(cell.sampling.extraction.tissue.name if cell.sampling else ''),
                        'Sampling Event': smart_text(cell.sampling.name if cell.sampling else ''),
                        'Group Color': str(color_map[cell]).replace('(', '[').replace(')', ']').replace(',', ''),
                        'Sampling Comment': smart_text(cell.sampling.comment if cell.sampling else ''),
                        'FACS Marker': smart_text(facs.marker.name if facs else ''),
                        'Cell Type': smart_text(cell.composition.name),
                        'Plate': smart_text(loc.plate.name),
                        'Well': smart_text(loc.well),
                        'Plate Location': smart_text(loc.plate.platestorage_set.all()[0] if loc.plate.platestorage_set.all() else '')
                    }


def print_cells_table(partner_name, individual_name=None, cell_folder=None, palette_name='hls'):
    cell_data_file = '{}cell_data.csv'.format(cell_folder)
    with open(cell_data_file, 'w') as f:
        fieldnames = ['CellContent ID'
                      'Cell ID',
                      'Cell Name',
                      'Cell Type',
                      'Plate',
                      'Well',
                      'Plate Location',
                      'Cell Group',
                      'Group Color',
                      'Sampling Event',
                      'Sampling Comment',
                      'Organ',
                      'Tissue',
                      'Sample Name',
                      'Sample Comment',
                      'Individual Name',
                      'Individual Comment',
                      'Gender',
                      'Sequencing File Name',
                      ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for cell_values in user_cells_table_values(partner_name, individual_name, cell_folder, palette_name=palette_name):
            writer.writerow(cell_values)
