import numpy as np
from scipy import stats
from math import sqrt, log10
from .hist import get_lims


def inflate_hist(hist, reads):
    pop = []
    for k in hist.keys():
        for i in range(int(hist[k]*reads)):
            pop.append(k)
    return pop


def pop_dist_sub(hist1, hist2):
    """
    Calculate the distance between two populations in the form of histograms
    Uses sum of deltas between columns hight divided by the number of columns
    """
    li, ri = get_lims(hist1, hist2)
    deltas = [abs(hist2[x]-hist1[x]) for x in range(li, ri)]
    score = sum(deltas)/(ri-li)
    return score


def pop_dist_subpeaks(hist1, hist2):
    """
    Calculate the distance between two populations in the form of histograms
    Uses sum of deltas beween columns hight divided by the number of columns
    """
    li, ri = get_lims(hist1, hist2)
    deltas = [abs(hist2[x]-hist1[x]) for x in range(li, ri)]
    score = sum(deltas) / (ri - li)
    return score


def pop_dist_ks_2samp(hist1, hist2, reads, sample_depth):
    """
    Calculate the distance between two populations in the form of histograms
    Uses 
    """
    pop1 = inflate_hist(hist1, reads)
    pop2 = inflate_hist(hist2, sample_depth)
    d, p = stats.ks_2samp(pop1, pop2)
    return 1-p


def alt_ks_2samp(hist1, hist2, reads, sample_depth):
    """
    Calculate the distance between two populations in the form of histograms
    Uses
    """
    pop1 = inflate_hist(hist1, reads)
    hist2.normalize()
    pop2 = np.random.choice(hist2._hist.keys(), sample_depth, p=hist2._hist.values())
    d, p = stats.ks_2samp(pop1, pop2)
    return 1-p


def alt2_ks_2samp(hist1, hist2, reads):
    """
    Calculate the distance between two populations in the form of histograms
    Uses
    """
    pop1 = inflate_hist(hist1, reads)
    hist2.normalize()
    pop2 = stats.rv_discrete(name='custm', values=(hist2._hist.keys(), hist2._hist.values()))
    d, p = stats.kstest(pop1, pop2.cdf)
    return 1-p


def pop_dist_corr(hist1, hist2):
    """
    Calculates the General Correlation Coefficient
    """
    li, ri = get_lims(hist1, hist2)
    dot_product = 0
    h1_sum_of_squares = 0
    h2_sum_of_squares = 0
    for x in range(li, ri):
        dot_product += hist1[x]*hist2[x]
        h1_sum_of_squares += hist1[x]**2
        h2_sum_of_squares += hist2[x]**2
    try:
        p = dot_product/sqrt(h1_sum_of_squares*h2_sum_of_squares)
    except ZeroDivisionError:
        print(type(hist1), type(hist2))
        print(hist1)
        print(hist2)
        raise
    return 1-p


def pop_dist_corr_numpy(hist1, hist2):
    """
    
    """
    li, ri = get_lims(hist1, hist2)
    pmat = np.corrcoef([hist1[x] for x in range(li, ri)], [hist2[x] for x in range(li, ri)])
    return 1-pmat[0][1]


def pop_dist_chisq(hist1, hist2):
    """
    
    """
    li, ri = get_lims(hist1, hist2)
    running_sum = 0
    for x in range(li, ri):
        running_sum += float((hist1[x]-hist2[x])**2)/(hist1[x]+hist2[x]) if hist1[x]+hist2[x]>0 else 0
    return 2*running_sum


def pop_dist_klp(hist1, hist2):
    """
    
    """
    li, ri = get_lims(hist1, hist2)
    running_sum = 0
    for x in range(li, ri):
        running_sum += hist1[x]*log10(hist1[x]/float(hist2[x]))
    return running_sum


def pop_dist_kl(hist1, hist2):
    return (pop_dist_klp(hist1, hist2) + pop_dist_klp(hist2, hist1))/float(2)


def prob(hist_sample, hist_distribution):
    li, ri = get_lims(hist_sample, hist_distribution)

    def zero_case_log(input_number):
        if input_number >= 0.001:
            return log10(input_number)
        return -3
    return -sum([zero_case_log(hist_distribution[bin])*hist_sample[bin] for bin in range(li, ri)])


def pop_dist(hist1, hist2, method='sub', reads=50, sample_depth=10000):
    """
    Calculate the distance between two populations in the form of histograms
    Method is given as a parameter 
    """
    if method == 'sub':
        return pop_dist_sub(hist1, hist2)
    if method == 'sp':
        return pop_dist_subpeaks(hist1, hist2)
    if method == 'ks':
        return pop_dist_ks_2samp(hist1, hist2, reads, sample_depth)
    if method == 'aks':
        return alt_ks_2samp(hist1, hist2, reads, sample_depth)
    if method == 'a2ks':
        return alt2_ks_2samp(hist1, hist2, reads)
    if method == 'cor':
        return pop_dist_corr(hist1, hist2)
    if method == 'con':
        return pop_dist_corr_numpy(hist1, hist2)
    if method == 'chi':
        return pop_dist_chisq(hist1, hist2)
    if method == 'kl':
        return pop_dist_kl(hist1, hist2)
    if method == 'pr':
        return prob(hist1, hist2)
    print('unknown method')
    raise