from collections import Counter


def shift_population(data, shift):
    """
    Shift population in the form of a histogram (Counter)
    """
    shifted_data = Counter()
    for i in list(data.keys()):
        shifted_data.update({i+shift: data[i]})
    return shifted_data


def shift_raw_population(data, shift):
    """
    Shift population in the form of raw values i.e. [14,14,15,15,15,16]
    """
    shifted_data = []
    for i in data:
        shifted_data.append(i+shift)
    return shifted_data


def normalize(counter_hist):
    """
    Normalize histogram (Counter)
    """
    m = float(sum(counter_hist.values()))
    for k in list(counter_hist.keys()):
        counter_hist[k] = counter_hist[k]/ m
    return counter_hist
