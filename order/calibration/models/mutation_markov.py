
import numpy

from .markov import FixedStepMarkovModel


class hashable_poly1d(numpy.poly1d):
    def __hash__(self):
        return hash((tuple(self.coeffs), self.order, self.variable))


def get_x_from_list(l, step):
    assert len(l) % step == 0
    for i in range(0, len(l), step):
        res = []
        for j in range(step):
            res.append(l[i+j])
        yield tuple(res)


class MutationMarkov(FixedStepMarkovModel):

    SQUEEZE = True

    def _steps(self, x):
        assert len(x) % 2 == 0
        # up_params = x[:len(x)//2]
        # dw_params = x[len(x)//2:]
        # ups = [hashable_poly1d([a, b, c]) for a, b, c in get_x_from_list(up_params, 3)]
        # dws = [hashable_poly1d([a, b, c]) for a, b, c in get_x_from_list(dw_params, 3)]
        up_params = x[:len(x)//4]
        dw_params = x[len(x)//4:]
        ups = [hashable_poly1d([a, b]) for a, b in get_x_from_list(up_params, 2)]
        dws = [hashable_poly1d([a, b]) for a, b in get_x_from_list(dw_params, 2)]
        # ups = [Hashable_exp(a, b) for a, b in get_x_from_list(up_params, 2)]
        # dws = [Hashable_exp(a, b) for a, b in get_x_from_list(dw_params, 2)]
        
        def fs(n):
            return 1-sum([fu(n) for fu in ups])-sum([fd(n) for fd in dws])
        
        d = {0: fs}
        d.update({i+1: fu for i, fu in enumerate(ups)})
        d.update({-i-1: fd for i, fd in enumerate(dws)})
        
        return d
