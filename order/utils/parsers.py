import csv
import io
import gzip

def parse_spaces_hist(hist, header):
    f = io.StringIO('{}\r\n{}'.format(header.strip(), hist.strip()))
    reader = csv.DictReader(f, delimiter=' ')
    row = next(reader)
    return {int(key) : int(row[key]) for key in list(row.keys())}


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


def uncalled_inputs(input_file, calling=None, reads_threshold=0):
    for loc, cell, row_hist in parse_input_file(input_file):
        if calling and calling[loc][cell]:
            continue
        if sum(row_hist.values()) < reads_threshold:
            continue
        yield loc, cell, row_hist