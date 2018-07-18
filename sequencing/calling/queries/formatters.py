from utils.user_cells_report import get_cells_color_map
from django.utils.encoding import smart_text
from sampling.models import FACS


def sr_labeler(sr):
    return 'ID{}'.format(sr.id)


def ms_labeler(ms):
    return 'LOC_{}'.format(ms.id)


def textify_keys_in_mutations_dict(mutation_dict, sr_label_func, ms_label_func):
    d = dict()
    for sr in mutation_dict:
        d[sr_label_func(sr)] = dict()
        for ms in mutation_dict[sr]:
            d[sr_label_func(sr)][ms_label_func(ms)] = mutation_dict[sr][ms]
    return d


def facs_or_none(sampling):
    if sampling:
        try:
            facs = sampling.facs
        except FACS.DoesNotExist:
            facs = None
    else:
        facs = None
    return facs


def get_cells_data_dict(srs, group_of_cell=lambda cell: '_'.join(cell.name.split('_')[:-2]), palette_name='hls', complet_color_map_override=None):
    if complet_color_map_override is None:
        groups = {
            g: i for i, g in enumerate(
                {
                    group_of_cell(cell) for cell in [
                        sr.barcoded_content.subclass.amplified_content.cell for sr in srs
                    ]
                 }
            )
        }
        cell_groups = {cell: groups[group_of_cell(cell)] for cell in
                       [sr.barcoded_content.subclass.amplified_content.cell for sr in srs]}
        color_map = get_cells_color_map(cell_groups, palette_name=palette_name)
    else:
        color_map = complet_color_map_override
    d = dict()
    for sr in srs:
        cell_cont = sr.barcoded_content.subclass.amplified_content
        cell = cell_cont.cell
        facs = facs_or_none(cell.sampling)
        for loc in cell_cont.physical_locations.exclude(plate__name__contains='AAR'):
            d[sr] = {
                'Content_ID': smart_text(sr.barcoded_content_id),
                'Cell_Group': group_of_cell(cell),
                'Cell_Group_short': group_of_cell(cell),
                'Cell_ID': smart_text(cell.pk),
                'Cell_Name': smart_text(cell.name),
                'Cell_Type': smart_text(cell.composition.name),
                'CellContent_ID': smart_text(cell_cont.pk),
                'Extraction_Event': smart_text(
                    cell.sampling.extraction.extraction_event.name if cell.sampling and cell.sampling.extraction and cell.sampling.extraction.extraction_event else ''),
                'FACS_Marker': smart_text(facs.marker.name if facs else ''),
                'Gender': smart_text(cell.individual.sex),
                'Group_Color': str(color_map[cell]).replace('(', '[').replace(')', ']').replace(',', '').replace(' ',
                                                                                                                 '_'),
                'Individual_Name': smart_text(cell.individual.name),
                'Organ': smart_text(cell.sampling.extraction.organ.name if cell.sampling else ''),
                'Plate': smart_text(loc.plate.name),
                'Plate_Location': smart_text(
                    loc.plate.platestorage_set.all()[0] if loc.plate.platestorage_set.all() else ''),
                'Sample_Name': smart_text(cell.sampling.extraction.name if cell.sampling else ''),
                'Sample_Reads_ID': smart_text(sr.id),
                'Sampling_Event': smart_text(cell.sampling.name if cell.sampling else ''),
                'Tissue': smart_text(cell.sampling.extraction.tissue.name if cell.sampling else ''),
                'Well': smart_text(loc.well),
            }
    return d
