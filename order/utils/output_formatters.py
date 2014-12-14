import csv
import gzip
from frogress import bar
from cloud.serialization.cloudpickle import loads, dumps

from parsers import parse_spaces_hist
from order.preprocessing import generate_sim_hists
from order.calling import call_hist


def get_or_create_call(hist, loc, cell, calling, sim_hists=None, min_cycles=0, max_cycles=100, method='cor', shift_margines=3, normalize=True, nsamples=None, trunc=False, cut_peak=False, trim_extremes=False):
    if calling[loc][cell]:
        return calling[loc][cell]
    if calling and sim_hists:
        calling[loc][cell] = call_hist(hist, sim_hists, min_cycles=min_cycles, max_cycles=max_cycles, method=method, shift_margines=shift_margines, normalize=normalize, nsamples=nsamples, trunc=trunc, cut_peak=cut_peak, trim_extremes=trim_extremes)
        return calling[loc][cell]
    print 'ERROR: No existing calling and not enough parameters to call from scratch'
    raise


def load_or_create_simulations_file(simulationsfile, max_cycles, normalize=True, truncate=False, cutpeak=False, trim_extremes=False):
    try:
        f = open(simulationsfile,'rb').read()
        print 'loading existing simulations'
        sim_hists = loads(f)
    except:
        print 'generating simulated histograms'
        sim_hists = generate_sim_hists(max_cycles, up=lambda d:0.00005*d**2 - 0.0009*d + 0.0036 , dw=lambda d:0.00009*d**2 - 0.00003*d - 0.0013, normalize=normalize, trunc=truncate, cut_peak=cutpeak, trim_extremes=trim_extremes) #15o
        #to#sim_hists = generate_sim_hists(200, up=lambda d:0.00002*d**2 - 0.0003*d + 0.0014 , dw=lambda d:0.00007*d**2 + 0.0005*d - 0.0031, normalize=normalize, trunc=truncate, cut_peak=cutpeak, trim_extremes=trim_extremes)
        with open(simulationsfile,'wb') as f:
            f.write(dumps(sim_hists))
        print 'done generating simulated histograms'
    return sim_hists


def generate_calling_file(input_file, 
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


def generate_output_file(input_file, output_file, calling,reads_threshold=50, score_threshold=0.006, verbose=False):
    with gzip.open(output_file, 'wb') as out:
        owrtr = csv.writer(out, dialect='excel-tab')
        with gzip.open(input_file, 'rb') as f:
            rdr = csv.reader(f, dialect='excel-tab')
            header_row = rdr.next()
            if verbose:
                header_row.extend(['shift', 'cycle', 'score', 'median', 'reads'])
            else:
                header_row.append('shift')
            owrtr.writerow(header_row)
            for row in bar(rdr):
                row_hist = parse_spaces_hist(row[2], header_row[2])
                loc = row[0]
                cell = row[1]
                if sum(row_hist.values()) < reads_threshold:
                    row.append('[]')
                else:
                    vc = get_or_create_call(row_hist, loc, cell, calling)
                    if vc['score'] > score_threshold:
                        row.append('[]')
                    else:
                        if verbose:
                            row.extend([vc['shift'], vc['cycle'], vc['score'], vc['median'], vc['reads']])
                        else:
                            row.append(vc['shift'])
                owrtr.writerow(row)