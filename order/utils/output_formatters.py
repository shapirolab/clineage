import csv
import gzip
from frogress import bar
from cloud.serialization.cloudpickle import loads, dumps

from parsers import parse_spaces_hist
from order.preprocessing import generate_sim_hists_of_up_to_k_alleles
from order.calling import call_multi_hist


def get_or_create_call(hist,
                       loc,
                       cell,
                       calling,
                       sim_hists,
                       **kwargs):
    """
        max_alleles=2,
        max_distance_from_median=30,

        min_cycles=0,
        max_cycles=100,
        method='cor',
        shift_margins=15,
        max_ms_length=60
        normalize=True,
        nsamples=None,
        truncate=False,
        cut_peak=False,
        trim_extremes=False
    :param hist:
    :param loc:
    :param cell:
    :param calling:
    :param sim_hists:
    :return:
    """
    if calling[loc][cell]:
        return calling[loc][cell]
    if calling and sim_hists:
        calling[loc][cell] = call_multi_hist(hist,
                                             sim_hists,
                                             **kwargs)
        return calling[loc][cell]
    print 'ERROR: No existing calling and not enough parameters to call from scratch'
    raise


def load_or_create_simulations_file(simulationsfile, **kwargs):
    """
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
        f = open(simulationsfile,'rb').read()
        print 'loading existing simulations'
        sim_hists = loads(f)
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
                          reads_threshold=50,
                          **kwargs):
    """
        max_alleles=2,
        max_distance_from_median=30,

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
    with gzip.open(input_file, 'rb') as f:
        rdr = csv.reader(f, dialect='excel-tab')
        header_row = rdr.next()
        for row in bar(rdr):
            row_hist = parse_spaces_hist(row[2], header_row[2])
            loc = row[0]
            cell = row[1]
            if sum(row_hist.values()) < reads_threshold:
                continue
            vc = get_or_create_call(row_hist, loc, cell, calling, sim_hists, **kwargs)
    return calling


def generate_output_file(input_file,
                         output_file,
                         calling,
                         sim_hists,
                         reads_threshold=50,
                         score_threshold=0.006,
                         verbose=False,
                         **kwargs):
    """
          shift_margins=3
          nsamples=None
          method='cor',
          score_threshold=0.006,
          min_cycles=0,
          max_cycles=80,
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
    print kwargs
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
                if sum(row_hist.values()) < reads_threshold or not calling[loc][cell]:
                    row.append('[]')
                else:
                    vc = get_or_create_call(row_hist, loc, cell, calling, sim_hists, score_threshold=score_threshold, **kwargs)
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