import os
import csv
from django.utils.encoding import smart_str
import seaborn as sns
from linapp.models import User, Cell, Individual


def query_partner_individuals(partner_name, individual_name=None):
    partner = User.objects.get(username__contains=partner_name)
    if individual_name:
        individuals = Individual.objects.filter(partner=partner).filter(name__contains=individual_name)
    else:
        individuals = Individual.objects.filter(partner=partner)
    return partner, individuals


def to_hex(n):
    return hex(int(n*255))[2:].upper()


def user_report_string(partner_name, individual_name=None):
    cellrow = 1
    color_map = get_cells_color_map(get_cells_grouping(partner_name, individual_name))
    report_string = '<!DOCTYPE html>\r\n<html>\r\n<body>\r\n<pre>\r\n'
    partner, individuals = query_partner_individuals(partner_name, individual_name)
    report_string += 'Collaborator: {}\r\n'.format(partner.username)
    for individual in individuals:
        report_string += '\tIndividual: {}\r\n'.format(individual.name)
        if not individual.extractionevent_set.all() or \
                not individual.extractionevent_set.all()[0].extraction_set.all() or \
                individual.extractionevent_set.all()[0].extraction_set.all()[0].samplingevent_set.all():
            report_string += '\t\t\tCells: {}\r\n'.format(individual.cell_set.count())
        report_string += '\tCollaborator table: {}\r\n'.format(None) #TODO:get this
        report_string += '\tDatabase table: {}\r\n'.format(None) #TODO:get this
        for ee in individual.extractionevent_set.all():
            for e in ee.extraction_set.all():
                report_string += '\t\tSample: {}\r\n'.format(e.name)
                report_string += '\t\tExtraction: {}\r\n'.format(ee.name)
                report_string += '\t\tExtraction date: {}\r\n'.format(ee.date)
                for se in e.samplingevent_set.all():
                    report_string += '\t\t\tCells separation: {}\r\n'.format(se.name)
                    report_string += '\t\t\tCells separation date: {}\r\n'.format(se.date)
                    report_string += '\t\t\tCells separation details: {}\r\n'.format(se.comment)
                    report_string += '\t\t\t<mark style="background-color:#{}{}{};">Cells: {}-{}</mark>\r\n'.format(to_hex(color_map[se.cell_set.all()[0]][0]),to_hex(color_map[se.cell_set.all()[0]][1]),to_hex(color_map[se.cell_set.all()[0]][2]),cellrow, cellrow+se.cell_set.count())
                    cellrow += se.cell_set.count()
    report_string += '</pre>\r\n</body>\r\n</html>\r\n'
    return report_string


def get_cell_files_from_folder(cell_folder):
    folder_cells = {}
    for file_name in os.listdir(cell_folder):
        if ('hist_' in file_name) or (not 'hist' in file_name):
            continue
        folder_cells[Cell.objects.get(id=int(file_name.split('-')[0]))] = file_name
    return folder_cells


def get_cells_filenames_in_folder(cells, cell_folder):
    folder_cells = get_cell_files_from_folder(cell_folder)
    for cell in cells:
        if cell in folder_cells:
            yield cell, folder_cells[cell]
        else:
            yield cell, None


def get_cells_grouping(partner_name, individual_name=None, current_group=0):
    partner, individuals = query_partner_individuals(partner_name, individual_name)
    cell_groups = {}
    for individual in individuals:
        if not individual.extractionevent_set.all() or \
                not individual.extractionevent_set.all()[0].extraction_set.all() or \
                not individual.extractionevent_set.all()[0].extraction_set.all()[0].samplingevent_set.all():
            for cell in individual.cell_set.all():
                cell_groups[cell] = current_group
            current_group += 1
            continue
        for ee in individual.extractionevent_set.all():
            for e in ee.extraction_set.all():
                for se in e.samplingevent_set.all():
                    for cell in se.cell_set.all():
                        cell_groups[cell] = current_group
                    current_group += 1
    return cell_groups


def get_cells_color_map(cell_groups):
    if len(set(cell_groups.values())) == 1:
        return {cell: (0,0,0) for cell in cell_groups}
    palette = sns.color_palette('hls', len(set(cell_groups.values())))
    return {cell: palette[cell_groups[cell]] for cell in cell_groups}


def user_cells_table_values(partner_name, individual_name=None, cell_folder=None):
    partner, individuals = query_partner_individuals(partner_name, individual_name)
    color_map = get_cells_color_map(get_cells_grouping(partner_name, individual_name))
    for individual in individuals:
        for cell, file_name in get_cells_filenames_in_folder(individual.cell_set.all(), cell_folder):
            for cell_cont in cell.cellcontent_set.all():
                assert cell_cont.physical_locations.exclude(plate__name__contains='AAR').count() <= 1
                for loc in cell_cont.physical_locations.exclude(plate__name__contains='AAR'):
                    yield {
                        'Cell ID': smart_str(cell.pk),
                        'Sequencing File Name': smart_str(file_name),
                        'Cell Name': smart_str(cell.name),
                        'Individual Name': smart_str(cell.individual.name),
                        'Individual Comment': smart_str(cell.individual.comment),
                        'Gender': smart_str(cell.individual.sex),
                        'Sample Name': smart_str(cell.sampling.extraction.name if cell.sampling else ''),
                        'Sample Comment': smart_str(cell.sampling.extraction.comment if cell.sampling else ''),
                        'Organ': smart_str(cell.sampling.extraction.organ.name if cell.sampling else ''),
                        'Tissue': smart_str(cell.sampling.extraction.tissue.name if cell.sampling else ''),
                        'Sampling Event': smart_str(cell.sampling.name if cell.sampling else ''),
                        'Group Color': color_map[cell],
                        'Sampling Comment': smart_str(cell.sampling.comment if cell.sampling else ''),
                        'Cell Type': smart_str(cell.composition.name),
                        'Plate': smart_str(loc.plate.name),
                        'Well': smart_str(loc.well)
                    }


def print_cells_table(partner_name, individual_name=None, cell_folder=None):
    cell_data_file = '{}cell_data.csv'.format(cell_folder)
    with open(cell_data_file, 'w') as f:
        fieldnames = ['Cell ID',
                      'Cell Name',
                      'Cell Type',
                      'Plate',
                      'Well',
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
        for cell_values in user_cells_table_values(partner_name, individual_name, cell_folder):
            writer.writerow(cell_values)