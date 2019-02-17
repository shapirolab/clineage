import glob
import os
import math
import sys
import decimal
import numpy as np
from collections import Counter


def get_lims(hist1, hist2):
    if not hist1.keys():
        if not hist2.keys():
            return 0, 0
        return min(hist2.keys()), max(hist2.keys()) + 1
    if not hist2.keys():
        return min(hist1.keys()), max(hist1.keys()) + 1
    li = min(min(hist1.keys()), min(hist2.keys()))
    ri = max(max(hist1.keys()), max(hist2.keys())) + 1
    if ri == li:
        ri = li + 1
    return int(li), int(ri)


def vnormalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2 == 0] = 1
    return a / np.expand_dims(l2, axis)


class Histogram(object):
    def __init__(self, h, normalize=False, nsamples=None, truncate=False, cut_peak=False, trim_extremes=False,
                 **kwargs):
        if isinstance(h, list):
            self._hist = Counter({i + 3: x for i, x in enumerate(h)})
        if isinstance(h, dict):
            self._hist = Counter(h)
        if isinstance(h, Histogram):
            self._hist = h._hist.copy()
            self.nsamples = h.nsamples
        else:
            if nsamples is not None:
                self.nsamples = nsamples
            else:
                self.nsamples = sum(self.values())
        self._vec = np.repeat(float(0), 50)
        for x in self._hist:
            self._vec[x] = self._hist[x]
        if trim_extremes:
            self.trim_extremes()
        if normalize:
            self.normalize()
        if truncate:
            self.truncate()
        if normalize:
            self.normalize()
        if cut_peak:
            self.cut_peak()
        if normalize:
            self.normalize()
        # self.clean_zero_entries()  # apparantly production runs faster without this

    # Cleaning
    def clean_zero_entries(self):
        for key in self.keys():
            if self[key] == 0:
                del self._hist[key]

    def trim_extremes(self, p=0.5):
        "Cleans extreme values over given percentage"
        extr_items = int(self.nsamples * p) if self.nsamples * p > 1.0 else 0
        self._hist = Counter(sorted(self.sample)[extr_items:-extr_items])

    def truncate(self, p=.050):
        "Cleans noise below p"
        for k in self.keys():
            if self[k] < p:
                self.nsamples -= self[k] * self.nsamples
                self[k] = 0

    def cut_peak(self, n=1):
        "Cleans anything with n zen zeros between it and the maximum"
        keys = self.keys()
        max_key = max([(self[k], k) for k in keys])[1]
        max_ind = keys.index(max_key)
        zeros_left = n
        for k in keys[max_ind:]:
            if not zeros_left:
                self.nsamples -= self[k] * self.nsamples
                self[k] = 0
            else:
                if self[k]:
                    zeros_left = n
                else:
                    zeros_left -= 1
        zeros_left = n
        for k in keys[max_ind::-1]:
            if not zeros_left:
                self.nsamples -= self[k] * self.nsamples
                self[k] = 0
            else:
                if self[k]:
                    zeros_left = n
                else:
                    zeros_left -= 1

    # Setters / Getters
    def keys(self):
        return sorted(self._hist.keys())

    def values(self):
        return [self._hist[k] for k in self.keys()]

    def items(self):
        return self._hist.items()

    def pop(self, item):
        return self._hist.pop(item)

    def __getitem__(self, item):
        return self._hist[item]

    def __setitem__(self, item, value):
        self._hist[item] = value

    @property
    def sample(self):
        return [k for k in self.keys()
                for i in range(int(self.nsamples * self[k]))]

    def random_sample(self, k):
        return Histogram(Counter(np.random.choice(self.keys(), k, p=self.values())))

    # Operators
    def normalize(self):
        # self.sq_normalize()
        s = float(sum(self.values()))
        if not s:
            return
        for k in self.keys():
            self._hist[k] /= s

    def sq_normalize(self, axis=-1, order=2):
        tuples_list = self._hist.items()
        keys = [t[0] for t in tuples_list]
        values = [t[1] for t in tuples_list]
        nvalues = vnormalized(values)[0]
        self._hist = Counter({k: v for k, v in zip(keys, nvalues)})

    def __add__(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i + other: self[i] for i in self.keys()}, nsamples=self.nsamples)
        if isinstance(other, Histogram):
            self.normalize()
            other.normalize()
            return Histogram({k: self[k] + other[k] for k in range(*get_lims(self, other))}, normalize=True)
        raise TypeError()

    def asym_add(self, other):
        if isinstance(other, Histogram):
            return Histogram({k: self[k] + other[k] for k in range(*get_lims(self, other))})
        raise TypeError()

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i - other: self[i] for i in self.keys()}, nsamples=self.nsamples)
        if isinstance(other, Histogram):
            self.normalize()
            other.normalize()
            return Histogram({k: self[k] - other[k] for k in range(*get_lims(self, other))}, normalize=True)

        raise TypeError()

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i * other: self[i] for i in self.keys()}, nsamples=self.nsamples)
        raise TypeError()

    def __div__(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i / other: self[i] for i in self.keys()}, nsamples=self.nsamples)
        raise TypeError()

    def __pow__(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i ** other: self[i] for i in self.keys()}, nsamples=self.nsamples)
        raise TypeError()

    def yadd(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i: self[i] + other for i in self.keys()}, nsamples=self.nsamples)
        if isinstance(other, decimal.Decimal):
            return self.yadd(float(other))
        raise TypeError()

    def ymul(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i: self[i] * other for i in self.keys()}, nsamples=self.nsamples)
        if isinstance(other, decimal.Decimal):
            return self.ymul(float(other))
        raise TypeError()

    def ydiv(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i: self[i] / other for i in self.keys()}, nsamples=self.nsamples)
        raise TypeError()

    def ypow(self, other):
        if isinstance(other, (int, float)):
            return Histogram({i: self[i] ** other for i in self.keys()}, nsamples=self.nsamples)
        raise TypeError()

    # Statistical operators
    def mu(self):
        return sum(k * self[k] for k in self.keys())

    def sig(self):
        return math.sqrt(mu((self - mu(self)) ** 2))

    def skew(self):
        if not sig(self):
            return 0
        return mu(((self - mu(self)) / sig(self)) ** 3)

    def copy(self, normalize=False):
        return Histogram(self._hist, nsamples=self.nsamples, normalize=normalize)

    # Repr
    def __repr__(self):
        N = sum(self.values())
        if N < .1:
            return "<Empty on [%s]>" % (', '.join('%.2f' % k for k in self.keys()))
        return '\n'.join(('%.2f: %.2f' % (x, self[x])).ljust(20)[:20] + '|'
                         + '#' * int(50 * self[x] / N + .5) for x in self.keys() if self[x])


def mu(x):
    return x.mu()


def sig(x):
    return x.sig()


def skew(x):
    return x.skew()
