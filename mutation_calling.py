import argparse
from order.calling import load_or_create_calling
from order.utils.output_formatters import generate_calling_file, \
    generate_output_file, \
    load_or_create_simulations_file, \
    save_calling_file


if '__main__' == __name__:
    parser = argparse.ArgumentParser(description='Analyses hist-pairs file')
    parser.add_argument('-n', '--normalize', type=bool, dest='normalize', default=True, help='path for tab delimited histpath-dup_id file')
    parser.add_argument('-t', '--truncate', type=bool, dest='truncate', default=False, help='Cleans noise below p (default p=.050)')
    parser.add_argument('-c', '--cutpeak', type=bool, dest='cutpeak', default=False, help='Cleans anything with n zen zeros between it and the maximum (default n=1)')
    parser.add_argument('-e', '--trimextremes', type=bool, dest='trimextremes', default=False, help='ignores extreme reads')
    parser.add_argument('-o', '--outputfile', type=str, dest='output', default='output.tab.gz', help='path of output file')
    parser.add_argument('-i', '--inputfile', type=str, dest='input', default='input.taqb.gz', help='path of input file')
    parser.add_argument('-r', '--minreads', type=int, dest='minreads', default=10, help='minimal reads threshold')
    parser.add_argument('-mc', '--maxcycles', type=int, dest='maxcycles', default=50, help='maximum cycles to simulate')
    parser.add_argument('-nc', '--mincycles', type=int, dest='mincycles', default=0, help='minimum cycles to simulate')
    parser.add_argument('-s', '--maxscore', type=float, dest='maxscore', default=999.9, help='max (worst) score allowed')
    parser.add_argument('-m', '--method', type=str, dest='method', default='emd', help='method for histogram comparison')
    parser.add_argument('-v', '--verbose', type=bool, dest='verbose', default=False, help='prints additional calling information columns in the output table')
    parser.add_argument('-sf', '--simulationsfile', type=str, dest='simulationsfile', default='sim_hists.pickle', help='path of preprocessed simulations file')
    parser.add_argument('-cf', '--callingfile', type=str, dest='callingfile', default='calling.pickle', help='path of preprocessed calling file')
    parser.add_argument('-sm', '--shiftmargins', type=int, dest='shiftmargins', default=15, help='number of attempted shifts to either side of the histogram median')
    parser.add_argument('-bd', '--simulationmethod', type=str, dest='simulationmethod', default='bin', help='method of MS histogram simulation')
    parser.add_argument('-w', '--workers', type=int, dest='workers', default=1, help='Number of working processes')
    #TODO: document those:
    parser.add_argument('-ma', '--maxalleles', type=int, dest='max_alleles', default=2, help='number of attempted shifts to either side of the histogram median')
    parser.add_argument('-sd', '--sampledepth', type=int, dest='sample_depth', default=10000, help='number of attempted shifts to either side of the histogram median')
    parser.add_argument('-ml', '--maxmslength', type=int, dest='max_ms_length', default=60, help='number of attempted shifts to either side of the histogram median')
    parser.add_argument('-md', '--mediandistance', type=int, dest='max_distance_from_median', default=10, help='number of attempted shifts to either side of the histogram median')

    args = parser.parse_args()
    normalize = args.normalize
    truncate = args.truncate
    cutpeak = args.cutpeak
    trim_extremes = args.trimextremes
    output_file = args.output
    input_file = args.input
    reads_threshold = args.minreads
    score_threshold = args.maxscore
    cycles_threshold = args.maxcycles
    min_cycles = args.mincycles
    shift_margins = args.shiftmargins
    max_alleles = args.max_alleles
    sample_depth = args.sample_depth
    max_ms_length = args.max_ms_length
    max_distance_from_median = args.max_distance_from_median
    meth = args.method
    verbose = args.verbose
    sim_method = args.simulationmethod
    SIMULATED_HISTS_PATH = args.simulationsfile
    SIGNALS_CALLING_PATH = args.callingfile
    workers = args.workers

    sim_hists = load_or_create_simulations_file(SIMULATED_HISTS_PATH,
                                                method=sim_method,
                                                max_cycles=cycles_threshold,
                                                up=lambda d: 0.00005*d**2 - 0.0009*d + 0.0036,
                                                dw=lambda d: 0.00009*d**2 - 0.00003*d - 0.0013,
                                                max_ms_length=max_ms_length,
                                                sample_depth=sample_depth,
                                                max_alleles=max_alleles,
                                                normalize=normalize,
                                                truncate=truncate,
                                                cut_peak=cutpeak,
                                                trim_extremes=trim_extremes)

    calling = load_or_create_calling(SIGNALS_CALLING_PATH)
    calling = generate_calling_file(input_file,
                                    sim_hists,
                                    calling,
                                    method=meth,
                                    reads_threshold=reads_threshold,
                                    workers=workers,
                                    score_threshold=score_threshold,
                                    min_cycles=min_cycles,
                                    max_cycles=cycles_threshold,
                                    nsamples=None,
                                    normalize=normalize,
                                    truncate=truncate,
                                    cut_peak=cutpeak,
                                    trim_extremes=trim_extremes,
                                    shift_margins=shift_margins,
                                    max_distance_from_median=max_distance_from_median,
                                    max_alleles=max_alleles,
                                    max_ms_length=max_ms_length
                                    )
    save_calling_file(calling, SIGNALS_CALLING_PATH)
    generate_output_file(input_file,
                         output_file,
                         calling,
                         reads_threshold=reads_threshold,
                         score_threshold=score_threshold,
                         verbose=verbose)