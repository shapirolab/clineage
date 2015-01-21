from scipy import optimize
import itertools
from itertools import repeat as rep
from multiprocessing import Pool
from time import sleep, time

import numpy
from order.hist import Histogram
from experimental_data.synthetic_data import get_hists_pairs, get_transposed_hists_pairs
from order.calibration.score import distance_from_model
from order.preprocessing import generate_hist


def doWork(input_tuple):
    alg, hist_pairs, cycles_tup = input_tuple
    cycles_26, cycles_47 = cycles_tup
    def mes(up, dw):
        score = 0.0
        for syn_len, syn_hist_26, syn_hist_47 in hist_pairs:
            msmodel = (up, dw, 'bin')
            score += distance_from_model(syn_len, syn_hist_26, cycles_26, msmodel, distance_measure='con')
            score += distance_from_model(syn_len, syn_hist_47, cycles_47, msmodel, distance_measure='con')
        return score
    def nmes(x):
        assert len(x) % 2 == 0
        up_params = list(x)[:len(x)/2]
        dw_params = list(x)[len(x)/2:]
        up = numpy.poly1d(up_params)
        dw = numpy.poly1d(dw_params)
        return mes(up, dw)
    
    # print nmes([0.0, 0.003, 0.0, 0.022])
    # minimizer_kwargs = dict(method="L-BFGS-B", bounds=[(0, 0.1), (-0.1, 0.1), (-0.5, 0.5), (0, 0.1), (-0.1, 0.1), (-0.5, 0.5)], options={'eps' : 1e-3, 'disp' : True})
    # minimizer_kwargs = dict(method="Nelder-Mead", bounds=[(0, 0.1), (-0.1, 0.1), (-0.1, 0.1), (-0.5, 0.5), (0, 0.1), (-0.1, 0.1), (-0.1, 0.1), (-0.5, 0.5)], options={'eps' : 1e-3, 'disp' : True})
    minimizer_kwargs = dict(method="Nelder-Mead", bounds=[(0, 0.1), (-0.1, 0.1), (-0.5, 0.5), (0, 0.1), (-0.1, 0.1), (-0.5, 0.5)], options={'eps' : 1e-3, 'disp' : True})
    # minimizer_kwargs = dict(method="TNC", bounds=[(0, 0.1), (-0.1, 0.1), (-0.5, 0.5), (0, 0.1), (-0.1, 0.1), (-0.5, 0.5)], options={'eps' : 1e-3, 'disp' : True})
    #res = optimize.basinhopping(nmes, [0.00005,  -0.0009, 0.0036, 0.00009, -0.00003, - .0013], minimizer_kwargs=minimizer_kwargs, niter=100)
    res = optimize.basinhopping(nmes, [4.83043475e-05, -8.84566705e-04, 3.63957747e-03, 9.75361458e-05, -2.80872372e-05, -1.24677989e-03], minimizer_kwargs=minimizer_kwargs, niter=2)
    return res


def inputs_generator():
    hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab')
    alg = 'con'
    for c1, c2 in itertools.product(range(11, 12, 1), range(38, 39, 1)):
        if c2 > c1:
            yield alg, hist_pairs, (c1, c2)


if __name__ == '__main__':
    # pool = Pool(processes=23)
    # rs = pool.map_async(doWork, inputs_generator(), chunksize=1)
    print len(list(inputs_generator()))
    rs = map(doWork, inputs_generator())
    print rs