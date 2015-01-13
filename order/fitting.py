from hist_dist import pop_dist
#import numpy as np
#from cv2 import *


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
