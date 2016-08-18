
from numpy import zeros
from numpy.linalg import matrix_power

from .base import ModelParams
from order.preprocessing import generate_mat_hist


class FixedStepMarkovModel(ModelParams):

    SQUEEZE = False

    def __init__(self, n=50, hist_kwargs=None):
        self.n = n
        self.hist_kwargs = hist_kwargs or {}
        self.hist_kwargs.update(dict(n=n))

    def _steps(self, x):
        raise NotImplementedError()

    def get_for_x(self, x):
        return self.Model(
            step_funcs=self._steps(x),
            n=self.n,
            squeeze=self.SQUEEZE,
            hist_kwargs=self.hist_kwargs,
        )

    class Model(object):

        def __init__(self, step_funcs, n, squeeze, hist_kwargs):
            mat = zeros((n, n))
            for step, func in step_funcs.items():
                for i in range(n):
                    val = func(i)
                    if 0 <= i + step < n:
                        mat[i, i + step] = val
                    elif squeeze:
                        ind = min(n-1, max(0, i+step))
                        mat[i, ind] += val
            self.mat = mat
            self.hist_kwargs = hist_kwargs

        def get_for_cycles(self, cycles):
            return self.CyclesModel(
                mat=self.mat,
                cycles=cycles,
                hist_kwargs=self.hist_kwargs,
            )

        class CyclesModel(object):

            def __init__(self, mat, cycles, hist_kwargs):
                self.mat = matrix_power(mat, cycles)
                self.hist_kwargs = hist_kwargs

            def get_hist_for_length(self, length):
                return generate_mat_hist(length, self.mat, **self.hist_kwargs)
