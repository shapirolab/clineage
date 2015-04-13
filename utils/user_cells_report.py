import os
import csv
from django.utils.encoding import smart_str
import seaborn as sns
from collections import defaultdict
from linapp.models import User, Cell, Individual, UserReport


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
        if not individual.extractionevent_set.all() or \
                not individual.extractionevent_set.all()[0].extraction_set.all() or \
                not individual.extractionevent_set.all()[0].extraction_set.all()[0].samplingevent_set.all():
            report_dict[partner.username][individual.name]['cells_list'] = 1 #temporary workaround TODO:revise
            for cls in sorted(list(set([cell.classification for cell in individual.cell_set.all()]))):
                report_dict[partner.username][individual.name][str(cls)]['Cells_color'] = hex_to_rgb(color_map, individual.cell_set.filter(classification=cls)[0])
                report_dict[partner.username][individual.name][str(cls)]['Cells_pos'] = [cellrow+1, cellrow+individual.cell_set.filter(classification=cls).count()]
                cellrow += individual.cell_set.filter(classification=cls).count()
        report_dict[partner.username][individual.name]['Collaborator_table'] = None
        if individuals.count() > 1:
            report_dict[partner.username][individual.name]['Database_table'] = 's:/LINEAGE/Hiseq/NSR2/fastq_human/Output/{}_cell_data.csv'.format(partner.username)
        else:
            report_dict[partner.username][individual.name]['Database_table'] = 's:/LINEAGE/Hiseq/NSR2/fastq_human/Output/{}_{}_cell_data.csv'.format(partner.username, individual.name)
        for ee in sorted(list(individual.extractionevent_set.all()), key=lambda i: i.name):
            for e in sorted(list(ee.extraction_set.all()), key=lambda i: i.name):
                report_dict[partner.username][individual.name][ee.name][e.name]['Extraction_date'] = ee.date
                for se in sorted(list(e.samplingevent_set.all()), key=lambda i: i.name):
                    for cls in sorted(list(set([cell.classification for cell in se.cell_set.all()]))):
                        if se.cell_set.filter(classification=cls):
                            report_dict[partner.username][individual.name][ee.name][e.name][se.name][str(cls)]['Cells_separation_date'] = se.date
                            report_dict[partner.username][individual.name][ee.name][e.name][se.name][str(cls)]['Cells_separation_details'] = se.comment
                            report_dict[partner.username][individual.name][ee.name][e.name][se.name][str(cls)]['Cells_color'] = hex_to_rgb(color_map, se.cell_set.filter(classification=cls)[0])
                            report_dict[partner.username][individual.name][ee.name][e.name][se.name][str(cls)]['Cells_pos'] = [cellrow+1, cellrow+se.cell_set.filter(classification=cls).count()]
                            cellrow += se.cell_set.filter(classification=cls).count()
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
    if not cell_groups:
        cell_groups = get_cells_grouping(partner_name, individual_name=individual_name)
    return cell_groups.keys()


def sorted_cells(individual):
    cell_list = []
    for ee in sorted(list(individual.extractionevent_set.all()), key=lambda ee: ee.name):
        for e in sorted(list(ee.extraction_set.all()), key=lambda e: e.name):
            for se in sorted(list(e.samplingevent_set.all()), key=lambda se: se.name):
                for cls in sorted(list(set([cell.classification for cell in se.cell_set.all()]))):
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
            for cls in sorted(list(set([cell.classification for cell in individual.cell_set.all()]))):
                for cell in individual.cell_set.filter(classification=cls):
                    cell_groups[cell] = current_group
                if individual.cell_set.filter(classification=cls):
                    current_group += 1
            continue
        for ee in sorted(list(individual.extractionevent_set.all()), key=lambda ee: ee.name):
            for e in sorted(list(ee.extraction_set.all()), key=lambda e: e.name):
                for se in sorted(list(e.samplingevent_set.all()), key=lambda se: se.name):
                    for cls in sorted(list(set([cell.classification for cell in se.cell_set.all()]))):
                        for cell in se.cell_set.filter(classification=cls):
                            cell_groups[cell] = current_group
                        if se.cell_set.filter(classification=cls):
                            current_group += 1
    return cell_groups


def get_cells_color_map(cell_groups, palette_name='hls'):
    if len(set(cell_groups.values())) == 1:
        return {cell: (0, 0.5, 0.5) for cell in cell_groups}
    palette = sns.color_palette(palette_name, len(set(cell_groups.values())))
    return {cell: palette[cell_groups[cell]] for cell in cell_groups}


def user_cells_table_values(partner_name, individual_name=None, cell_folder=None, palette_name='hls'):
    partner, individuals = query_partner_individuals(partner_name, individual_name)
    color_map = get_cells_color_map(get_cells_grouping(partner_name, individual_name), palette_name)
    for individual in sorted(list(individuals), key=lambda i: i.name):
        for cell, file_name in get_cells_filenames_in_folder(sorted_cells(individual), cell_folder):
            for cell_cont in cell.cellcontent_set.all():
                assert cell_cont.physical_locations.exclude(plate__name__contains='AAR').count() <= 1
                for loc in cell_cont.physical_locations.exclude(plate__name__contains='AAR'):
                    yield {
                        'Cell ID': smart_str(cell.pk),
                        'Sequencing File Name': smart_str(file_name),
                        'Cell Name': smart_str(cell.name),
                        'Cell Group': smart_str(cell.classification),
                        'Individual Name': smart_str(cell.individual.name),
                        'Individual Comment': smart_str(cell.individual.comment),
                        'Gender': smart_str(cell.individual.sex),
                        'Sample Name': smart_str(cell.sampling.extraction.name if cell.sampling else ''),
                        'Sample Comment': smart_str(cell.sampling.extraction.comment if cell.sampling else ''),
                        'Organ': smart_str(cell.sampling.extraction.organ.name if cell.sampling else ''),
                        'Tissue': smart_str(cell.sampling.extraction.tissue.name if cell.sampling else ''),
                        'Sampling Event': smart_str(cell.sampling.name if cell.sampling else ''),
                        'Group Color': list(color_map[cell]),
                        'Sampling Comment': smart_str(cell.sampling.comment if cell.sampling else ''),
                        'Cell Type': smart_str(cell.composition.name),
                        'Plate': smart_str(loc.plate.name),
                        'Well': smart_str(loc.well)
                    }


def print_cells_table(partner_name, individual_name=None, cell_folder=None, palette_name='hls'):
    cell_data_file = '{}cell_data.csv'.format(cell_folder)
    with open(cell_data_file, 'w') as f:
        fieldnames = ['Cell ID',
                      'Cell Name',
                      'Cell Type',
                      'Plate',
                      'Well',
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