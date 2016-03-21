
import csv
from io import StringIO
from utils.wells import xy_index_2_index, index2str
import re
from sampling.models import SampleComposition
from django.utils.encoding import smart_bytes

def plate_parser(plate_string, rows_number=8, columns_number=12):
    """
    Convert an copy-from-excel plate representation (tab delimited values) to a list of (well_index, well value) tuples.
    *Sorted according to plate-global well index
    Example:
    	1	2	3	4	5	6	7	8	9	10	11	12
    A	1	2		4		6	7	8	9	10	11	12
    B	1	2	3	4	5	6	7	8	9	10	11	12
    C	1		3	4	5	6	7	8	9	10	11	12
    D	1	2	3	4	5	6	7	8	9	10	11	12
    E	1	2	3	4	5	6	7	8	9	10	11	12
    F	1	2	3	4	5	6	7	8	9	10	11	12
    G	1	2	3	4	5	6	7	8	9	10	11	12
    H	1	2	b1	b2	b3	b4	b5	b6	b7	b8

                        ||
                        ||
                        \/

[(1, '1'), (2, '1'), (3, '1'), (4, '1'), (5, '1'), (6, '1'), (7, '1'), (8, '1'), (9, '2'), (10, '2'), (11, ''), (12, '2'), (13, '2'), (14, '2'), (15, '2'), (16, '2'), (17, ''), (18, '3'), (19, '3'), (20, '3'), (21, '3'), (22, '3'), (23, '3'), (24, 'b1'), (25, '4'), (26, '4'), (27, '4'), (28, '4'), (29, '4'), (30, '4'), (31, '4'), (32, 'b2'), (33, ''), (34, '5'), (35, '5'), (36, '5'), (37, '5'), (38, '5'), (39, '5'), (40, 'b3'), (41, '6'), (42, '6'), (43, '6'), (44, '6'), (45, '6'), (46, '6'), (47, '6'), (48, 'b4'),
 (49, '7'), (50, '7'), (51, '7'), (52, '7'), (53, '7'), (54, '7'), (55, '7'), (56, 'b5'), (57, '8'), (58, '8'), (59, '8'), (60, '8'), (61, '8'), (62, '8'), (63, '8'), (64, 'b6'), (65, '9'), (66, '9'), (67, '9'), (68, '9'), (69, '9'), (70, '9'), (71, '9'), (72, 'b7'), (73, '10'), (74, '10'), (75, '10'), (76, '10'), (77, '10'), (78, '10'), (79, '10'), (80, 'b8'), (81, '11'), (82, '11'), (83, '11'), (84, '11'), (85, '11'), (86, '11'), (87, '11'), (88, ''), (89, '12'), (90, '12'), (91, '12'), (92, '12'), (93, '12'), (9
4, '12'), (95, '12'), (96, '')]
    """
    plate_string.replace('\x00', '')
    plate_string = smart_bytes(plate_string)
    plate_string_io = StringIO(plate_string)
    dialect = csv.Sniffer().sniff(plate_string)  # TODO: default to 'excel-tab' dialect
    has_header = csv.Sniffer().has_header(plate_string)
    reader = list(csv.reader(plate_string_io, dialect))
    indexed_values = []

    assert rows_number <= len(reader) <= rows_number + 1
    # try:
    for row_index, row in enumerate(reader):
        assert columns_number <= len(row) <= columns_number + 1
        for col_index, col in enumerate(row):
            if row_index == 0 and col_index == 0 and not row[col_index]:  # handle labels as top row
                has_header = True
                break
            if has_header and col_index == 0 and re.match('[a-z]+$', row[col_index], re.IGNORECASE):
                continue
            if has_header:  # in the header case the indices are already one-based
                well_index = xy_index_2_index(row_index, col_index, dimentions=(rows_number, columns_number))
            else:  # in the other (no header) case the indices are zero-based so we add 1
                well_index = xy_index_2_index(row_index+1, col_index+1, dimentions=(rows_number, columns_number))
            indexed_values.append((well_index, row[col_index]))
    return sorted(indexed_values, key=lambda tuple_value: tuple_value[0])
    # except csv.Error as e:
    #     return 'line %d: %s' % (reader.line_num, e)  # TODO: return appropriate error data.

def well_value_to_composition(well_value):
    if re.match('[0-9]+$', well_value, re.IGNORECASE):
        if int(well_value) == 1:
            return SampleComposition.objects.get(name='Single Cell')  # single_cell
        if int(well_value) > 1:
            return SampleComposition.objects.get(name='Bulk')  #bulk
    if well_value.lower() in ['b', 'bulk']:
        return SampleComposition.objects.get(name='Bulk')  #bulk
    return None

def test():
    # TODO: write some tests

    print('Test passed successfully')


if __name__ == '__main__':
    test()