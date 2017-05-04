
import operator


def highest_peak(hist):
    """
    Simple highest peak calling.
    The proportion of the highest peak from the rest of the histogram is provided as confidence
    Args:
        hist:

    Returns:

    """
    key, value = max(hist.items(), key=operator.itemgetter(1))
    confidence = value/sum(hist.values())
    return key, confidence


def highest_peak_diff_from_second(hist):
    """
    return the key of the
    Args:
        hist:

    Returns:

    """
    if len(hist) == 1:
        key, value = hist.popitem()
        return key, 1.0
    sorted_items = sorted(hist.items(), key=lambda x: x[1], reverse=True)
    (key, max_value), (_, second_value) = sorted_items[:2]
    return key, (max_value-second_value)/sum(hist.values())



