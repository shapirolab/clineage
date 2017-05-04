
import numpy
from pandas.core.algorithms import isin

from .markov import FixedStepMarkovModel


def get_x_from_list(l, step):
    assert len(l) % step == 0
    for i in range(0, len(l), step):
        res = []
        for j in range(step):
            res.append(l[i+j])
        yield tuple(res)


class PropagationMarkov(FixedStepMarkovModel):

    SQUEEZE = False

    def __init__(self, p_degree, degrees_dict, *args, **kwargs):
        self.p_degree = p_degree
        self.degrees_dict = degrees_dict
        super().__init__(*args, **kwargs)

    def _calculate_polynomes(self, x):
        assert len(x) % 2 == 0
        p = numpy.poly1d(x[:self.p_degree])
        d = dict()
        for step in self.degrees_dict:
            start, end = self.degrees_dict[step]
            d[step] = numpy.poly1d(x[start:end])
        return p, d

    def _calculate_steps(self, p, d):
        assert 0 not in d.keys()
        fs = 1 - sum(list(d.values()))
        if isinstance(p, numpy.poly1d):
            pd = {0: 1 + (fs * p)}
            pd.update({k: f * p for k, f in d.items()})
        elif isinstance(p, list):
            pd = {0: lambda x: 1 + (fs(x) * p[x])}
            pd.update({k: lambda x, f=f: f(x) * p[x] for k, f in list(d.items())})
            # We add f=f to avoid late binding
        return pd

    def _steps(self, x):
        p, d = self._calculate_polynomes(x)
        pd = self._calculate_steps(p, d)
        return pd


class PropagationMarkovConstantP(PropagationMarkov):
    def __init__(self, p, *args, **kwargs):
        self.p = p
        super().__init__(*args, **kwargs)

    def _calculate_polynomes(self, x):
        assert len(x) % 2 == 0
        d = dict()
        for step in self.degrees_dict:
            start, end = self.degrees_dict[step]
            d[step] = numpy.poly1d(x[start:end])
        return self.p, d


# class PropagationMarkovStepP(PropagationMarkov):
#     def __init__(self, p, *args, **kwargs):
#         self.p = p
#         super().__init__(*args, **kwargs)
#
#     class Model(object):
#
#         def __init__(self, step_funcs, n, squeeze, hist_kwargs):
#             p_x, pry = PropagationMarkovStepP.p
#             mat = numpy.zeros((n, n), dtype='float128')
#             for step, func in step_funcs.items():
#                 for i in range(n):
#                     val = func(i)
#                     if 0 <= i + step < n:
#                         mat[i, i + step] = val
#                     elif squeeze:
#                         ind = min(n-1, max(0, i+step))
#                         mat[i, ind] += val
#             self.mat = mat
#             self.hist_kwargs = hist_kwargs
