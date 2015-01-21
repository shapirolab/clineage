from scipy import optimize
import itertools
from itertools import repeat as rep
from multiprocessing import Pool
from time import sleep, time

import numpy
from experimental_data.synthetic_data import get_hists_pairs, get_transposed_hists_pairs
from order.calibration.score import distance_from_model, distance_from_model_across_lengths


def doWork(input_tuple):
    alg, optimizer_method, hist_pairs, cycles_tup = input_tuple
    cycles_26, cycles_47 = cycles_tup
    def nmes(x):
        assert len(x) % 2 == 0
        up_params = list(x)[:len(x)/2]
        dw_params = list(x)[len(x)/2:]
        up = numpy.poly1d(up_params)
        dw = numpy.poly1d(dw_params)
        msmodel = (up, dw, optimizer_method)
        return distance_from_model_across_lengths(msmodel, hist_pairs, cycles_tup, distance_measure=alg)

    minimizer_kwargs = dict(method=optimizer_method, bounds=[(0, 0.1), (-0.1, 0.1), (-0.5, 0.5), (0, 0.1), (-0.1, 0.1), (-0.5, 0.5)], options={'eps' : 1e-3, 'disp' : True})
    res = optimize.basinhopping(nmes, [0.00005,  -0.0009, 0.0036, 0.00009, -0.00003, -0.0013], minimizer_kwargs=minimizer_kwargs, niter=100)
    return res


def inputs_generator():
    hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab')
    alg = 'con'
    optimizer_method = "L-BFGS-B"
    # optimizer_method = "Nelder-Mead"
    # optimizer_method = "TNC"
    for c1, c2 in itertools.product(range(11, 12, 1), range(38, 39, 1)):
        if c2 > c1:
            yield alg, optimizer_method, hist_pairs, (c1, c2)


if __name__ == '__main__':
    # pool = Pool(processes=23)
    # rs = pool.map_async(doWork, inputs_generator(), chunksize=1)
    print len(list(inputs_generator()))
    rs = map(doWork, inputs_generator())
    print rs