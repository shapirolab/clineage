
from order.calibration.score import distance_from_model_across_lengths
from scipy import optimize


class Struct(object):

    def __init__(self, **kwargs):
        for n in self.__slots__:
            setattr(self, n, kwargs[n])


class Optimizer(Struct):

    __slots__ = [
        "minimizer_method",
        "dist_metric",
        "iterations",
        "disp",
        # "minimizer_kwargs",
        # "basinhopping_kwargs",
    ]

    def optimize_across_lengths(self, params, steps_func, reference, **kwargs):
        def nmes(x):
            steps = steps_func(x)
            return distance_from_model_across_lengths(steps, reference,
                self.dist_metric, **kwargs)
        return optimize.basinhopping(
            func=nmes,
            x0=params.initial_guess,
            niter=self.iterations,
            T=params.temperature,
            stepsize=params.stepsize,
            minimizer_kwargs=dict(
                method=self.minimizer_method,
                bounds=params.bounds,
                options=dict(
                    disp=self.disp,
                    **params.optimizer_options
                ),
            ),
            disp=self.disp,
        )


class Params(Struct):

    __slots__ = [
        "bounds",
        "initial_guess",
        "stepsize",
        "temperature",
        "optimizer_options",
        # "minimizer_kwargs",
    ]


class Reference(Struct):

    __slots__ = [
        "syn_hist_list",
        "cycles_tup",
    ]
