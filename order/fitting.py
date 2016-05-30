from .hist_dist import pop_dist


def match_cycles(hist, sim_hists, method='cor', reads=50, min_cycles=0, max_cycles=100, **kwargs):
    s = 9999
    c = 0
    best_sim_hist = None
    for cycles in range(min_cycles, max_cycles):
        sim_hist = sim_hists[cycles]
        score = pop_dist(hist, sim_hist, method=method, reads=reads)
        if score < s:
            s = score
            c = cycles
            best_sim_hist = sim_hist
    return c, s, best_sim_hist
