import csv
import gzip
from frogress import bar
from cloud.serialization.cloudpickle import loads, dumps

from parsers import parse_spaces_hist
from order.preprocessing import generate_sim_hists_of_up_to_k_alleles
from order.calling import generate_hist_calls

def load_or_create_simulations_file(simulationsfile, **kwargs):
    """
        method='bin',
        max_cycles=90,
        up=lambda d:0.00005*d**2 - 0.0009*d + 0.0036,
        dw=lambda d:0.00009*d**2 - 0.00003*d - 0.0013,
        normalize=True,
        truncate=False,
        cut_peak=False,
        trim_extremes=False,
        max_ms_length=60,
        sample_depth=10000,
        max_alleles = 2

    :param simulationsfile:
    :param kwargs:
    :return:
    """
    try:
        f = open(simulationsfile, 'rb').read()
        print 'loading existing simulations'
        sim_hists = loads(f)
        print 'done loading existing simulations'
    except:
        print 'generating simulated histograms'
        sim_hists = generate_sim_hists_of_up_to_k_alleles(**kwargs)
        with open(simulationsfile, 'wb') as f:
            f.write(dumps(sim_hists))
        print 'done generating simulated histograms'
    return sim_hists


def generate_calling_file(input_file,
                          sim_hists,
                          calling,
                          **kwargs):
    """
        workers=1,
        max_alleles=2,
        max_distance_from_median=30,
        reads_threshold=50
        shift_margins=3
        nsamples=None
        method='cor',
        score_threshold=0.006,
        min_cycles=0,
        max_cycles=80,
        max_ms_length=60
        normalize=True,
        truncate=False,
        cut_peak=False,
        trim_extremes=False):
    :param input_file:
    :param sim_hists:
    :param calling:
    :param kwargs:
    :return:
    """
    for loc, cell, row_hist, res in bar(generate_hist_calls(input_file, sim_hists, calling, **kwargs)):
        calling[loc][cell] = res
    return calling


def generate_output_file(input_file,
                         output_file,
                         calling,
                         reads_threshold=50,
                         score_threshold=0.006,
                         verbose=False):
    """
    :param input_file:
    :param sim_hists:
    :param calling:
    :return:
    """
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
                loc = row[0]
                cell = row[1]
                if not calling[loc][cell] or calling[loc][cell]['reads'] < reads_threshold:
                    row.append('[]')
                else:
                    vc = calling[loc][cell]
                    if vc['score'] > score_threshold:
                        row.append('[]')
                    else:
                        if verbose:
                            row.extend([str(vc['shifts']), vc['cycle'], vc['score'], vc['median'], vc['reads']])
                        else:
                            row.append(vc['shifts'])
                owrtr.writerow(row)


def save_calling_file(calling, callingfile):
    try:
        f = open(callingfile, 'rb').read()
    except:
        with open(callingfile, 'wb') as f:
            f.write(dumps(calling))