import pacal
from collections import defaultdict
from order.hist import Histogram


def sim(g, up, dw):
    x = pacal.BinomialDistr(g, up)
    y = pacal.BinomialDistr(g, dw)
    return x-y


def update_probs(cycle, up, dw, probs):
    if cycle <= 0:
        return probs
    new_probs = defaultdict(float)
    for key in probs.keys():
        up_p = up(key)
        dw_p = dw(key)
        new_probs[key] += (1-up_p-dw_p)*probs[key]
        new_probs[key-1] += dw_p*probs[key]
        new_probs[key+1] += up_p*probs[key]
    return update_probs(cycle-1, up, dw, Histogram(new_probs, normalize=False, nsamples=1, truncate=False, cut_peak=False, trim_extremes=False) + probs)
        
    
def dyn_prob(cycle, d, up, dw, normalize=False, nsamples=1, truncate=False, cut_peak=False, trim_extremes=False):
    hist = Histogram({d:1.0}, normalize=normalize, nsamples=1, truncate=truncate, cut_peak=cut_peak, trim_extremes=trim_extremes)
    return update_probs(cycle, up, dw, hist)