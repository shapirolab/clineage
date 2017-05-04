import numpy as np
from .hist import get_lims


def pop_dist_corr_numpy(hist1, hist2):
    """
    1-correlation
    """
    li, ri = get_lims(hist1, hist2)
    pmat = np.corrcoef([hist1[x] for x in range(li, ri)], [hist2[x] for x in range(li, ri)])
    return 1-pmat[0][1]
