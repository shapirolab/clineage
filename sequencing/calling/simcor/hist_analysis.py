import numpy as np
import peakutils


def remove_points_close_to_top(hist, hs, number_of_points, distance=1):
    """
    Removes points on the histogram if they are within the given range, and returns the histogram without them
    Searches if there's a point in the given range, and removes it,
    :param hist: histogram
    :param hs: sorted histogram
    :param distance: the distance to remove
    :return: new histogram
    """
    for steps_from_point in range(1, distance):
        if hist[hs[number_of_points][0] + steps_from_point] > 0:
            hist.pop(hs[number_of_points][0] + steps_from_point)
        if hist[hs[number_of_points][0] - steps_from_point] > 0:
            hist.pop(hs[number_of_points][0] - steps_from_point)
    return hist


def get_far_apart_highest_peaks(hist, allele_number=1, minimal_distance_between_peaks=1, min_prop=0):
    """
    Identify the k highest peaks that satisfy minimal distance
    Args:
        hist: histogram
        allele_number: allele_number
        minimal_distance_between_peaks: minimum distance between the allele
        min_prop: peaks with proportion lower than this are skipped 
    """
    histogram = hist.copy()
    histogram.normalize()
    for allele in range(allele_number):
        hs = sorted([(k, v) for k, v in histogram.items() if v > min_prop], key=lambda hkey: hkey[1], reverse=True)
        if allele >= len(hs):
            break
        histogram = remove_points_close_to_top(histogram, hs, allele, minimal_distance_between_peaks)
    hs = sorted([(k, v) for k, v in histogram.items() if v > min_prop], key=lambda hkey: hkey[1], reverse=True)
    seeds = [x for x, y in hs[:allele_number]]
    return seeds


def better_get_far_apart_highest_peaks(hist, minimal_distance_between_peaks=1, min_prop=0.02):
    """
    Identify the k highest peaks that satisfy minimal distance
    Args:
        hist: histogram
        allele_number: allele_number
        minimal_distance_between_peaks: minimum distance between the allele
        min_prop: peaks with proportion lower than this are skipped 
    """
    histogram = hist.copy()
    histogram.normalize()
    x = list(range(min(histogram.keys())-5, max(histogram.keys()) + 6))  # add zero padding so that peakutils won't get stuck
    v = np.array([histogram[i] for i in x])
    indexes = peakutils.indexes(v, thres=min_prop/max(v), min_dist=minimal_distance_between_peaks)
    return [x[i] for i in indexes]


import numpy as np
import scipy.signal
def better_get_far_apart_highest_peaks_that_doesnt_hang(hist, minimal_distance_between_peaks=1):
    """
    Identify the k highest peaks that satisfy minimal distance
    Args:
        hist: histogram
        allele_number: allele_number
        minimal_distance_between_peaks: minimum distance between the allele
        min_prop: peaks with proportion lower than this are skipped 
    """
    histogram = hist.copy()
    histogram.normalize()
    x = list(range(min(histogram.keys()), max(histogram.keys()) + 1))
    v = np.array([histogram[i] for i in x])
    peak_widths = np.arange(1, 2*minimal_distance_between_peaks+1)
    peak_indices = scipy.signal.find_peaks_cwt(v, peak_widths)
    return [x[i] for i in peak_indices]
