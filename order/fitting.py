import numpy as np
from histogram_utils import normalize, shift_population
from collections import Counter, defaultdict
from hist_dist import pop_dist
#from cv2 import *


def inflate_hist(hist, reads):
    pop = []
    for k in hist.keys():
        for i in range(int(hist[k]*reads)):
            pop.append(k)
    return pop

# def convert_to_cv():
    # a = np.zeros((5,2))
    # for i in range(0,5):
        # a[i][1] = i+1
    
    # a[0][0] = 1
    # a[1][0] = 1
    # a[2][0] = 0
    # a[3][0] = 0
    # a[4][0] = 1

    # b = np.zeros((4,2))

    # for i in range(0,4):
        # b[i][1] = i+1

    # b[0][0] = 0
    # b[1][0] = 1
    # b[2][0] = 0
    # b[3][0] = 1

    # a64 = cv.fromarray(a)
    # a32 = cv.CreateMat(a64.rows, a64.cols, cv.CV_32FC1)
    # cv.Convert(a64, a32)
    # b64 = cv.fromarray(b)
    # b32 = cv.CreateMat(b64.rows, b64.cols, cv.CV_32FC1)
    # cv.Convert(b64, b32)

    # cv.CalcEMD2(a32,b32,cv.CV_DIST_L2)


def match_cycles(hist, sim_hists, ms_len_for_probabilities, method='cor', reads=50, min_cycles=0, max_cycles=100):
    s = 9999
    c = 0
    best_sim_hist = None
    for cycles in range(min_cycles, max_cycles):
        sim_hist = sim_hists[ms_len_for_probabilities][cycles]
        score = pop_dist(hist, sim_hist, method=method, reads=reads)
        if score < s:
            s = score
            c = cycles
            best_sim_hist = sim_hist
    return c, s, best_sim_hist


def fit(ms_hist, sim_hists, method='cor' ,reads=50, min_cycles=0, max_cycles=100):
    ml = np.median(ms_hist.sample)
    out = []
    for d in range(int(np.median(ms_hist.sample))-3, int(np.median(ms_hist.sample))+3):
        normalized_shifted_reads_hist = normalize(shift_population(Counter(ms_hist.sample), -d))  
        c, s, best_sim_hist = match_cycles(normalized_shifted_reads_hist, sim_hists, d, method=method ,reads=reads, min_cycles=min_cycles, max_cycles=max_cycles)
        out.append([d, c, s, best_sim_hist])
    locis = defaultdict(lambda: defaultdict(list))
    res = out[0]
    for r in out:
        if r[2] < res[2]:
            res = r
    return res
