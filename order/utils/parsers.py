import csv
import io
import gzip
from frogress import bar
from collections import defaultdict

def parse_spaces_hist(hist, header):
    f = io.StringIO('{}\r\n{}'.format(header.strip(), hist.strip()))
    reader = csv.DictReader(f, delimiter=' ')
    row = next(reader)
    return {int(key) : int(row[key]) for key in row.keys()}


def parse_input_file(input_file):
    """
    :param input_file:
    :param reads_threshold:
    :return:
    """
    with gzip.open(input_file, 'rb') as f:
        rdr = csv.reader(f, dialect='excel-tab')
        header_row = next(rdr)
        for row in rdr:
            row_hist = parse_spaces_hist(row[2], header_row[2])
            loc = row[0]
            cell = row[1]
            yield loc, cell, row_hist


def uncalled_inputs(input_file, calling, reads_threshold=0):
    for loc, cell, row_hist in parse_input_file(input_file):
        if calling[loc][cell]:
            continue
        if sum(row_hist.values()) < reads_threshold:
            continue
        yield loc, cell, row_hist

        
def parse_output_file(output_file):
    calling = defaultdict(lambda: defaultdict(dict))
    with gzip.open(output_file, 'rb') as f:
        rdr = csv.DictReader(f, dialect='excel-tab')
        for row in bar(rdr):
            if row['shift'] == '[]':
                continue
            calling[row['loc_name']][row['cell_path']] = {'hist': parse_spaces_hist(row['3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150'], '3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150'), 'shifts': eval(row['shift']), 'cycle':int(row['cycle']), 'score':float(row['score']), 'median':int(row['median']), 'reads':int(row['reads']), }
    return calling