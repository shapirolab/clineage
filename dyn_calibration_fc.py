import itertools
from experimental_data.synthetic_data import get_hists_pairs, get_transposed_hists_pairs
from order.calibration.score import optimize_across_lengths
import argparse
from pickle import dumps

def inputs_generator(cycles_range):
    """
    cycles_range=((i,j),(k,l))
    """
    c1, c2 = cycles_range
    hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab', lengths=[5, 10, 15, 20, 25, 30])
    #hist_pairs = get_transposed_hists_pairs(filename='experimental_data/hist_by_ms_len_as_0_sum.tab', lengths=[30])
    alg = 'con'
    optimizer_method = "L-BFGS-B"
    bounds = [(0.0, 0.0001), (-0.005, 0.005), (-0.05, 0.05), (0.0, 0.0001), (-0.005, 0.005), (-0.05, 0.05)]
    optimizer_options = {'eps': 1e-5, 'disp': False}
    # initial_guess = [0.00005,  -0.0009, 0.0036, 0.00009, -0.00003, -0.0013]
    initial_guess = [0.0,  0.0, 0.0, 0.0, 0.0, 0.0]
    iterations = 100
    sim = 'mat'
    # optimizer_method = "Nelder-Mead"
    # optimizer_method = "TNC"
    # optimizer_method = "L-BFGS-B"
    for c1, c2 in itertools.product(range(*c1), range(*c2)):
        if c2 > c1:
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
    print c1, c2
    for result in map(optimize_across_lengths, inputs_generator(cycles_range=((c1, c1+1), (c2, c2+1)))):
        results.append(result)
        alg, sim, optimizer_method, cycles_tup, ms_lengths, res = result
        print alg, sim, optimizer_method, cycles_tup, ms_lengths, res
    with open('out/{}_{}_{}.pickle'.format(unique_id, n, user_id),'wb') as f:
        f.write(dumps(results))
    