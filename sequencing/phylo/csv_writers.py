import csv


def print_mutation_dict_to_file(textual_mutation_dict, output_path):
    loci_columns = {msl for sr in textual_mutation_dict for msl in textual_mutation_dict[sr]}
    with open(output_path, 'w') as f:
        fieldnames = ['names'] + list(loci_columns)
        writer = csv.DictWriter(f, fieldnames=fieldnames, dialect='excel-tab')
        writer.writeheader()
        for sr in textual_mutation_dict:
            row_dict = {'names': sr}
            row_dict.update(textual_mutation_dict[sr])
            row_dict.update({loc: 'NaN' for loc in loci_columns - textual_mutation_dict[sr].keys()})
            writer.writerow(row_dict)


def write_cell_data_dict_to_file(cell_data_dict, output_path):
    with open(output_path, 'w') as f:
        fieldnames = ['Content_ID', 'Cell_Group', 'Cell_Group_short', 'Cell_ID', 'Cell_Name', 'Cell_Type',
                      'CellContent_ID', 'Extraction_Event', 'FACS_Marker', 'Gender', 'Group_Color', 'Individual_Name',
                      'Organ', 'Plate', 'Plate_Location', 'Sample_Name', 'Sample_Reads_ID', 'Sampling_Event', 'Tissue',
                      'Well']
        writer = csv.DictWriter(f, fieldnames=fieldnames, dialect='excel-tab')
        writer.writeheader()
        for sr in cell_data_dict:
            writer.writerow(cell_data_dict[sr])
