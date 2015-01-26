import itertools
import concurrent.futures

from experimental_data.synthetic_data import get_hists_pairs, get_transposed_hists_pairs
from order.calibration.score import optimize_across_lengths


def inputs_generator():
    hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab', lengths=[5, 10, 15, 20, 25, 30])
    #hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab', lengths=[30])
    alg = 'con'
    optimizer_method = "L-BFGS-B"
    bounds = [(0, 0.001), (-0.01, 0.01), (-0.05, 0.05), (0, 0.1), (-0.1, 0.1), (-0.05, 0.05)]
    optimizer_options = {'eps' : 1e-5, 'disp' : False}
    # initial_guess = [0.00005,  -0.0009, 0.0036, 0.00009, -0.00003, -0.0013]
    initial_guess = [0.0,  0.0, 0.0, 0.0, 0.0, 0.0]
    iterations=100
    sim='bon'
    # optimizer_method = "Nelder-Mead"
    # optimizer_method = "TNC"
    # optimizer_method = "L-BFGS-B"
    for c1, c2 in itertools.product(range(10, 13, 1), range(42, 46, 1)):
    # for c1, c2 in itertools.product(range(11, 12, 1), range(38, 39, 1)):
        if c2 > c1:
            yield alg, sim, optimizer_method, hist_pairs, (c1, c2), bounds, initial_guess, iterations, optimizer_options


if __name__ == '__main__':
    # pool = Pool(processes=23)
    # rs = pool.map_async(doWork, inputs_generator(), chunksize=1)
    # print len(list(inputs_generator()))
    # rs = map(optimize_across_lengths, inputs_generator())
    # print rs
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
        # for result in map(optimize_across_lengths, inputs_generator()):
        for result in executor.map(optimize_across_lengths, inputs_generator()):
            results.append(result)
            alg, sim, optimizer_method, cycles_tup, ms_lengths, res = result
            print alg, sim, optimizer_method, cycles_tup, ms_lengths, res