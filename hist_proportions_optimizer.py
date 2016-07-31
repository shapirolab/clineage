from order.preprocessing import generate_simulated_proportional_alleles_precalculated
from order.hist_dist import pop_dist
import numpy as np
from scipy import optimize

optimizer_method = "L-BFGS-B"
bounds = [(0.0, 1.0)]
optimizer_options = {'eps': 1e-3, 'disp': False}


def fit_proportions(input_tuple):
    hist, seeds, sim_alleles, cycles_tup, dist_method, = input_tuple
    measured_hist_shift = int(np.mean(list(seeds)))
    normalized_shifted_reads_hist = hist-measured_hist_shift

    def nmes(p):
        sim_hist = generate_simulated_proportional_alleles_precalculated(seeds,
                                                                         sim_alleles,
                                                                         cycles_tup,
                                                                         (p, 1 - p))
        score = pop_dist(normalized_shifted_reads_hist, sim_hist, method=dist_method)
        return score

    res = optimize.fminbound(nmes, 0.0, 1.0)
    return res, seeds, cycles_tup


