from scipy import optimize
import itertools
from itertools import repeat as rep
from multiprocessing import Pool
from time import sleep, time

import numpy
from experimental_data.synthetic_data import get_hists_pairs, get_transposed_hists_pairs
from order.calibration.score import distance_from_model, distance_from_model_across_lengths, optimize_across_lengths


def inputs_generator():
    hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab')
    alg = 'con'
    optimizer_method = "L-BFGS-B"
    bounds = [(0, 0.1), (-0.1, 0.1), (-0.5, 0.5), (0, 0.1), (-0.1, 0.1), (-0.5, 0.5)]
    optimizer_options = {'eps' : 1e-3, 'disp' : True}
    initial_guess = [0.00005,  -0.0009, 0.0036, 0.00009, -0.00003, -0.0013]
    iterations=100
    # optimizer_method = "Nelder-Mead"
    # optimizer_method = "TNC"
    for c1, c2 in itertools.product(range(11, 12, 1), range(38, 39, 1)):
        if c2 > c1:
            yield alg, optimizer_method, hist_pairs, (c1, c2), bounds, initial_guess, iterations, optimizer_options


if __name__ == '__main__':
    # pool = Pool(processes=23)
    # rs = pool.map_async(doWork, inputs_generator(), chunksize=1)
    print len(list(inputs_generator()))
    rs = map(optimize_across_lengths, inputs_generator())
    print rs