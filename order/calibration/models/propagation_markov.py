
import numpy

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

    def __init__(self, p_degree, degrees_dict, *kwargs):
        self.p_degree = p_degree
        self.degrees_dict = degrees_dict
        super().__init__(kwargs)

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
        fs = 1 - sum(d.values())
        pd = {0: 1 + (p * fs)}
        pd.update({k: p * f for k,f in d.items()})
        return pd

    def _steps(self, x):
        p, d = self._calculate_polynomes(x)
        pd = self._calculate_steps(p, d)
        return pd

