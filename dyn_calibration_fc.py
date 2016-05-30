import itertools
from experimental_data.synthetic_data import get_hists_pairs, get_transposed_hists_pairs
from order.calibration.score import optimize_across_lengths
import argparse
from pickle import dumps
import concurrent.futures

def inputs_generator(cycles_range):
    """
    cycles_range=((i,j),(k,l))
    """
    c1, c2 = cycles_range
    hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab', lengths=[5, 10, 15, 20, 25, 30])
    # hist_pairs = get_transposed_hists_pairs(filename='experimental_data/AG_hist_by_ms_len_as_0_sum.tab', lengths=[5, 10, 15, 20, 25])
    # hist_pairs = get_transposed_hists_pairs(filename='experimental_data/A_hist_by_ms_len_as_0_sum.tab', lengths=[5, 10, 19])
    #hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab', lengths=[30])
    alg = 'con'
    # optimizer_method = "L-BFGS-B"

    # bounds = [(0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1),
    #           (0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1),
    #           (0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1),
    #           (0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1),
    #           (0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1),
    #           (0.000, 0.0001), (-0.01, 0.01), (-0.1, 0.1)]
    # bounds = [(-0.005, 0.005), (-0.01, 0.01),
    #           (-0.005, 0.005), (-0.01, 0.01),
    #           (-0.005, 0.005), (-0.01, 0.01),
    #           (-0.005, 0.005), (-0.01, 0.01),
    #           (-0.005, 0.005), (-0.01, 0.01),
    #           (-0.005, 0.005), (-0.01, 0.01),
    #           (-0.005, 0.005), (-0.01, 0.01),
    #           (-0.005, 0.005), (-0.01, 0.01),]
    bounds = [(0.00, 0.005), (-0.1, 0.1),
              (0.00, 0.005), (-0.1, 0.1),
              (0.00, 0.005), (-0.1, 0.1),
              (0.00, 0.005), (-0.1, 0.1),]
    optimizer_options = {'eps': 1e-5, 'disp': True}
    # initial_guess = [  1.03705828e-04,   1.82528119e-04,
    #                  1.76578108e-05,   1.81881361e-04,
    #                  2.34019025e-05,   2.61905091e-04,
    #                  0.0,  0.0,
    #                  1.80806879e-03,   -1.00000000e-02,
    #                  2.77248056e-04,   -1.82022699e-03,
    #                  8.04505476e-05,   -9.21648191e-04,
    #                  7.41017931e-05,   -5.89879970e-04]
    # initial_guess = [2.35880857e-04,  -1.63107341e-03,
    #                 1.77516360e-03,  -9.68678589e-03,
    #                 -3.66578315e-05,  -1.93883174e-04,
    #                 0, 0,
    #                 ]
    initial_guess = [0, 0,
                     0, 0,
                     0, 0,
                     0, 0,]
    #                  0, 0,
    #                  0, 0,
    #                  0, 0,
    #                  0, 0,]
    # initial_guess = [0, 2.35880857e-04,  -1.63107341e-03,
    #                  0, 0, 0,
    #                  0, 0, 0,
    #                  0.0,  1.77516360e-03,  -9.68678589e-03,
    #                  1.45317052e-05,  -3.66578315e-05,  -1.93883174e-04,
    #                  0, 0, 0,]
    # initial_guess = [0, 0, 0,
    #                  0, 0, 0,
    #                  # 0, 0, 0,
    #                  # 0, 0, 0,
    #                  0, 0, 0,
    #                  0, 0, 0,]
    optimizer_method = "L-BFGS-B"
    iterations = 100
    sim = 'mat'
    # optimizer_methods = ["Nelder-Mead", "TNC", "Powell", "COBYLA", "SLSQP", "L-BFGS-B"]
    # optimizer_methods = ["COBYLA", "SLSQP", "L-BFGS-B"]
    # optimizer_method = "Nelder-Mead"
    # optimizer_method = "TNC"
    # optimizer_method = "Powell"
    # optimizer_method = "COBYLA"
    optimizer_methods = ["L-BFGS-B"]
    # optimizer_method = "L-BFGS-B"

    for c1, c2 in itertools.product(list(range(*c1)), list(range(*c2))):
        if c2 > c1:
            for optimizer_method in optimizer_methods:
                yield alg, sim, optimizer_method, hist_pairs, (c1, c2), bounds, initial_guess, iterations, optimizer_options


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyses hist-pairs file')
    parser.add_argument('-n', '--run_number', type=int, dest='run_number', default=1, help='run number from which the cycles range is deduced')
    parser.add_argument('-i', '--run_id', type=int, dest='run_id', default=1, help='run id')
    parser.add_argument('-u', '--user_id', type=str, dest='user_id', default='user', help='user id')
    parser.add_argument('-fn', '--c1_min', type=int, dest='c1_min', default=1, help='c1 cycles range')
    parser.add_argument('-sn', '--c2_min', type=int, dest='c2_min', default=1, help='c2 cycles range')
    parser.add_argument('-fm', '--c1_max', type=int, dest='c1_max', default=1, help='c1 cycles range')
    parser.add_argument('-sm', '--c2_max', type=int, dest='c2_max', default=1, help='c2 cycles range')
    args = parser.parse_args()
    n = args.run_number
    unique_id = args.run_id
    user_id = args.user_id
    c1_min = args.c1_min
    c2_min = args.c2_min
    c1_max = args.c1_max
    c2_max = args.c2_max
    
    assert n/(c1_max-c1_min) <= c2_max-c2_min
    c1 = c1_min + n % (c1_max-c1_min)
    c2 = c2_min + n / (c1_max-c1_min)
    results = []
    print(c1, c2)
    with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
        # for result in executor.map(optimize_across_lengths, inputs_generator(cycles_range=((c1, c1+1), (c2, c2+1)))):
        for result in map(optimize_across_lengths, inputs_generator(cycles_range=((c1, c1+1), (c2, c2+1)))):
            results.append(result)
            alg, sim, optimizer_method, cycles_tup, ms_lengths, res = result
            print(alg, sim, optimizer_method, cycles_tup, ms_lengths, res)
    with open('out/{}_{}_{}.pickle'.format(unique_id, n, user_id), 'wb') as f:
        f.write(dumps(results))