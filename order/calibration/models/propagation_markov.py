
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

    def _steps(self, x):
        assert len(x) % 2 == 0
        p = numpy.poly1d(x[:2]) 
        x = x[2:]
        # up_params = x[:len(x)//2]
        # dw_params = x[len(x)//2:]
        # ups = [hashable_poly1d([a, b, c]) for a, b, c in get_x_from_list(up_params, 3)]
        # dws = [hashable_poly1d([a, b, c]) for a, b, c in get_x_from_list(dw_params, 3)]
        up_params = x[:len(x)//4]
        dw_params = x[len(x)//4:]
        ups = [numpy.poly1d([a, b]) for a, b in get_x_from_list(up_params, 2)]
        dws = [numpy.poly1d([a, b]) for a, b in get_x_from_list(dw_params, 2)]
        # ups = [Hashable_exp(a, b) for a, b in get_x_from_list(up_params, 2)]
        # dws = [Hashable_exp(a, b) for a, b in get_x_from_list(dw_params, 2)]

        fs = 1 - sum(ups) - sum(dws)
        
        d = {0: 1+(p*fs)}
        d.update({i+1: p*fu for i, fu in enumerate(ups)})
        d.update({-i-1: p*fd for i, fd in enumerate(dws)})
        
        return d
