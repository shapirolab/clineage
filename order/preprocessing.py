from collections import defaultdict, Counter
from histogram_utils import normalize
from binomial_sim import sim
from frogress import bar as tqdm
from order.hist import Histogram

sample_depth = 10000

def generate_sim_hists(max_cycles, up=lambda x: 0.003, dw=lambda x: 0.022, normalize=True, nsamples=None, trunc=False, cut_peak=False, trim_extremes=False):
    sim_hists = defaultdict(dict)
    for d in tqdm(range(0,60)):
        for cycles in range(max_cycles):
            up_p = up(d)
            dw_p = dw(d)
            z = sim(cycles, up_p, dw_p)
            sim_hists[d][cycles] = Histogram(Counter(z.rand(sample_depth)), normalize=normalize, nsamples=sample_depth, trunc=trunc, cut_peak=cut_peak, trim_extremes=trim_extremes)
    return sim_hists