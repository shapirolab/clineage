import csv
import StringIO


def parse_spaces_hist(hist,header):
    f = StringIO.StringIO('{}\r\n{}'.format(header.strip(), hist.strip()))
    reader = csv.DictReader(f, delimiter=' ')
    row = reader.next()
    return {int(key) : int(row[key]) for key in row.keys()}

def parse_input_file(input_file, 
                      sim_hists, 
                      calling, 
                      method='cor', 
                      reads_threshold=50, 
                      score_threshold=0.006, 
                      min_cycles=0, 
                      max_cycles=80, 
                      normalize=True, 
                      truncate=False, 
                      cutpeak=False, 
                      trim_extremes=False):
    with gzip.open(input_file, 'rb') as f:
        rdr = csv.reader(f, dialect='excel-tab')
        header_row = rdr.next()
        for row in bar(rdr):
            row_hist = parse_spaces_hist(row[2], header_row[2])
            loc = row[0]
            cell = row[1]
            if sum(row_hist.values()) < reads_threshold:
                continue
            vc = get_or_create_call(row_hist, loc, cell, calling, sim_hists, min_cycles=min_cycles, max_cycles=max_cycles, method=method, shift_margines=3, normalize=True, nsamples=None, trunc=truncate, cut_peak=cutpeak, trim_extremes=trim_extremes)
    return calling
